// One three-body encounter: drawing the initial conditions, integrating, and
// classifying what is left afterwards.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <random>

#include "cluster.hpp"
#include "constants.hpp"
#include "three_body.hpp"
#include "vec3.hpp"

namespace bhs {

using Rng = std::mt19937_64;

// The standard fixes what the generator returns but leaves the distributions
// to each library, and the libraries differ: the same seed gave other speeds
// under GCC than under clang.  The distributions are therefore written out
// here, so that a seed means the same initial conditions everywhere.  They
// return what libc++ returns, with which the first results were made.

/// Uniform in [0, 1): one number of the generator.
inline double canonical(Rng& rng) { return static_cast<double>(rng()) / 18446744073709551616.0; }   // 2^64

inline double uniform(Rng& rng, double lo, double hi) { return (hi - lo) * canonical(rng) + lo; }

/// A whole number in [0, count), each as likely as the others: as many bits
/// of the generator as count - 1 has, drawn again if they are too large.
inline std::size_t uniformIndex(Rng& rng, std::size_t count) {
    if (count <= 1) return 0;
    int bits = 0;
    while (bits < 64 && ((count - 1) >> bits) != 0) ++bits;
    const std::uint64_t mask = bits == 64 ? ~std::uint64_t{0} : (std::uint64_t{1} << bits) - 1;
    std::uint64_t value;
    do {
        value = rng() & mask;
    } while (value >= count);
    return static_cast<std::size_t>(value);
}

/// Two independent Gaussians of unit dispersion (Marsaglia's polar method).
struct GaussianPair {
    double first, second;
};
inline GaussianPair gaussianPair(Rng& rng) {
    double u, v, s;
    do {
        u = uniform(rng, -1.0, 1.0);
        v = uniform(rng, -1.0, 1.0);
        s = u * u + v * v;
    } while (s > 1 || s == 0);
    const double factor = std::sqrt(-2 * std::log(s) / s);
    return {u * factor, v * factor};
}

/// Speed drawn from a Maxwellian of one-dimensional dispersion sigma: the
/// norm of three Gaussians.  (The fourth of the two pairs is not used.)
inline double maxwellSpeed(Rng& rng, double sigma) {
    const GaussianPair xy = gaussianPair(rng);
    const double x = xy.first * sigma, y = xy.second * sigma, z = gaussianPair(rng).first * sigma;
    return std::sqrt(x * x + y * y + z * z);
}

// ----------------------------------------------------------------------------
// Kepler orbits
// ----------------------------------------------------------------------------

/// Eccentric anomaly E solving E - e sin E = meanAnomaly on [0, 2 pi]
/// (Newton iteration safeguarded by bisection).
inline double solveKepler(double meanAnomaly, double e) {
    double lo = 0.0, hi = kTwoPi, x = meanAnomaly;
    for (int iter = 0; iter < 200; ++iter) {
        const double f = x - e * std::sin(x) - meanAnomaly;
        (f > 0 ? hi : lo) = x;
        double next = x - f / (1.0 - e * std::cos(x));
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        if (std::fabs(next - x) < 1e-15 * (1.0 + std::fabs(next))) return next;
        x = next;
    }
    return x;
}

/// How an orbit lies in space [rad].
struct Orientation {
    double node = 0.0, inclination = 0.0, argPericentre = 0.0;
};

inline Orientation randomOrientation(Rng& rng) {
    Orientation o;
    o.inclination = std::acos(uniform(rng, -1.0, 1.0));
    o.node = uniform(rng, 0.0, kTwoPi);
    o.argPericentre = uniform(rng, 0.0, kTwoPi);
    return o;
}

/// Position and velocity of one body relative to another.
struct RelativeState {
    Vec3 r, v;
};

/// A point of a Kepler orbit with specific angular momentum h and
/// eccentricity e, at distance r and true anomaly f from the focus.
inline RelativeState orbitPoint(double mu, double h, double e, double r, double trueAnomaly, const Orientation& o) {
    const double theta = o.argPericentre + trueAnomaly;   // argument of latitude
    const double cn = std::cos(o.node), sn = std::sin(o.node);
    const double ci = std::cos(o.inclination), si = std::sin(o.inclination);
    const double ct = std::cos(theta), st = std::sin(theta);
    const Vec3 radial{cn * ct - sn * st * ci, sn * ct + cn * st * ci, st * si};
    const Vec3 transverse{-cn * st - sn * ct * ci, -sn * st + cn * ct * ci, ct * si};
    return {r * radial, (mu / h * e * std::sin(trueAnomaly)) * radial + (h / r) * transverse};
}

/// A bound orbit (a, e) at the given mean anomaly.
inline RelativeState boundOrbit(double mu, double a, double e, double meanAnomaly, const Orientation& o) {
    const double eccAnomaly = solveKepler(meanAnomaly, e);
    const double trueAnomaly = std::atan2(std::sqrt(1.0 - e * e) * std::sin(eccAnomaly), std::cos(eccAnomaly) - e);
    return orbitPoint(mu, std::sqrt(mu * a * (1.0 - e * e)), e, a * (1.0 - e * std::cos(eccAnomaly)), trueAnomaly, o);
}

/// Specific angular momentum and eccentricity of an unbound orbit with speed
/// vInf at infinity and pericentre rp.
struct UnboundOrbit {
    double h, ecc;
};

inline UnboundOrbit unboundOrbit(double mu, double vInf, double rp) {
    const double h = rp * std::sqrt(vInf * vInf + 2.0 * mu / rp);
    return {h, std::sqrt(1.0 + vInf * vInf * h * h / (mu * mu))};
}

/// Time that orbit takes between its pericentre and `distance` from the
/// focus, on the way in or on the way out.
inline double timeFromPericentre(double mu, double vInf, double rp, double distance) {
    if (vInf < 1e-12) {   // parabola
        const double d = std::sqrt(distance / rp - 1.0);
        return std::sqrt(2.0 * rp * rp * rp / mu) * (d + d * d * d / 3.0);
    }
    // hyperbola: r = aHyp (e cosh H - 1)
    const double aHyp = mu / (vInf * vInf);
    const double ecc = unboundOrbit(mu, vInf, rp).ecc;
    const double H = std::acosh((distance / aHyp + 1.0) / ecc);
    return std::sqrt(aHyp * aHyp * aHyp / mu) * (ecc * std::sinh(H) - H);
}

/// An unbound orbit (speed vInf at infinity, pericentre rp) on its way in, at
/// `distance` from the focus, and the time it needs from there to pericentre.
/// The way from infinity to `distance` is Kepler's, so nothing is integrated
/// before this point.
struct Approach {
    RelativeState state;
    double timeToPericentre = 0.0;
};

inline Approach incomingOrbit(double mu, double vInf, double rp, double distance, const Orientation& o) {
    const UnboundOrbit orbit = unboundOrbit(mu, vInf, rp);
    const double trueAnomaly =
        -std::acos(std::clamp((orbit.h * orbit.h / (mu * distance) - 1.0) / orbit.ecc, -1.0, 1.0));
    return {orbitPoint(mu, orbit.h, orbit.ecc, distance, trueAnomaly, o), timeFromPericentre(mu, vInf, rp, distance)};
}

/// Relative orbit of two bodies ("PairOrbit" in the notebooks).
struct PairOrbit {
    double energy = 0.0;          // specific orbital energy; negative if bound
    double semiMajorAxis = 0.0;   // infinite if unbound
    double eccentricity = 0.0;
    Vec3 separation, relVelocity, angMomentum;
};

inline PairOrbit pairOrbit(const Vec3& ri, const Vec3& vi, double mi, const Vec3& rj, const Vec3& vj, double mj) {
    PairOrbit orbit;
    orbit.separation = rj - ri;
    orbit.relVelocity = vj - vi;
    const double mu = kG * (mi + mj);
    const double r = norm(orbit.separation);
    orbit.energy = dot(orbit.relVelocity, orbit.relVelocity) / 2 - mu / r;
    orbit.semiMajorAxis = orbit.energy < 0 ? -mu / (2 * orbit.energy) : std::numeric_limits<double>::infinity();
    orbit.angMomentum = cross(orbit.separation, orbit.relVelocity);
    const Vec3 eccVector = cross(orbit.relVelocity, orbit.angMomentum) / mu - orbit.separation / r;
    orbit.eccentricity = norm(eccVector);
    return orbit;
}

// ----------------------------------------------------------------------------
// Initial conditions
// ----------------------------------------------------------------------------

/// One of the three bodies.  Stars all have the same mass and radius, and
/// can be pulled apart by the tide of a black hole.
struct Body {
    double mass = 0.0;
    bool isStar = false;
};
inline Body blackHole(double mass) { return {mass, false}; }
inline Body star() { return {kStarMass, true}; }

using Bodies = std::array<Body, 3>;

/// Everything needed to integrate one encounter.  Times count from the start
/// of the integration.
struct EncounterSetup {
    Bodies bodies{};
    Masses masses{};
    State initialState{};
    double approachTime = 0.0;  // the incoming orbit reaches its pericentre
    double firstLook = 0.0;     // the outcome is looked at for the first time
    double duration = 0.0;      // the integration ends, if the encounter is not settled before
    BinaryScales binary;        // the binary before the encounter
    double vInfinity = 0.0;     // speed of the incoming body at infinity
    double pericentre = 0.0;    // pericentre distance of the incoming orbit
    double inclination = 0.0;   // of the incoming orbit to the x-y plane [rad]
};

inline EncounterSetup newSetup(const Bodies& bodies, double mBinary, double a, double e) {
    EncounterSetup setup;
    setup.bodies = bodies;
    setup.masses = {bodies[0].mass, bodies[1].mass, bodies[2].mass};
    setup.binary = {a, e, kTwoPi * std::sqrt(a * a * a / (kG * mBinary))};
    return setup;
}

/// Tidal radius of a binary (semi-major axis a) that meets a single body:
/// closer than this, the tide of the single body pulls the binary apart.
inline double binaryTidalRadius(double a, double mBinary, double mSingle) {
    return a * std::cbrt(mSingle / mBinary);
}

/// How far from a binary a single body is taken to move on a Kepler orbit
/// about the binary's centre of mass: kStartTidalRadii tidal radii, and as
/// many semi-major axes where the single body is the lighter one.  The
/// incoming body starts here.
inline double startDistance(double a, double mBinary, double mSingle) {
    return kStartTidalRadii * std::max(a, binaryTidalRadius(a, mBinary, mSingle));
}

/// The times of an encounter whose incoming orbit (about the mass mu / G)
/// starts `approachTime` before its pericentre at the distance `start`, and
/// which is given up `latest` after the pericentre.
inline void setTimes(EncounterSetup& setup, double mu, double start, double approachTime, double latest) {
    setup.approachTime = approachTime;
    setup.duration = approachTime + latest;
    setup.firstLook = std::min(
        setup.duration, approachTime + timeFromPericentre(mu, setup.vInfinity, setup.pericentre, kFirstLook * start));
}

/// Whether the node and argument of pericentre of the incoming orbit are
/// drawn at random or left at zero.
enum class ApproachAngles { Zero, Random };

/// A binary (bodies 0 and 1; orbit a, e) meets body 2, which arrives from a
/// population of velocity dispersion sigma with a pericentre drawn uniformly
/// in [kMinPericentre, kMaxPericentre] a.  The binary lies in the x-y plane
/// with random orientation and phase, its centre of mass at rest at the
/// origin; the incoming orbit is inclined isotropically.  (Code A, B, C.)
inline EncounterSetup drawBinarySingle(const Bodies& bodies, double a, double e, double sigma, ApproachAngles angles,
                                       Rng& rng) {
    const double m1 = bodies[0].mass, m2 = bodies[1].mass, m3 = bodies[2].mass;
    EncounterSetup setup = newSetup(bodies, m1 + m2, a, e);
    setup.vInfinity = maxwellSpeed(rng, sigma);
    setup.pericentre = uniform(rng, kMinPericentre * a, kMaxPericentre * a);

    setup.inclination = std::acos(uniform(rng, -1.0, 1.0));

    Orientation inPlane;
    inPlane.argPericentre = uniform(rng, 0.0, kTwoPi);
    inPlane.node = uniform(rng, 0.0, kTwoPi);
    const double meanAnomaly = uniform(rng, 0.0, kTwoPi);
    const RelativeState inner = boundOrbit(kG * (m1 + m2), a, e, meanAnomaly, inPlane);   // from body 0 to body 1
    const Vec3 r1 = (m1 / (m1 + m2)) * inner.r, v1 = (m1 / (m1 + m2)) * inner.v;
    setBody(setup.initialState, 0, r1 - inner.r, v1 - inner.v);
    setBody(setup.initialState, 1, r1, v1);

    Orientation approach;
    approach.inclination = setup.inclination;
    if (angles == ApproachAngles::Random) {
        approach.node = uniform(rng, 0.0, kTwoPi);
        approach.argPericentre = uniform(rng, 0.0, kTwoPi);
    }
    const double mu = kG * (m1 + m2 + m3);
    const double start = startDistance(a, m1 + m2, m3);
    const Approach outer = incomingOrbit(mu, setup.vInfinity, setup.pericentre, start, approach);
    setBody(setup.initialState, 2, outer.state.r, outer.state.v);

    // In, and at the latest until a body leaving as it came is kEndDistance a away.
    setTimes(setup, mu, start, outer.timeToPericentre,
             timeFromPericentre(mu, setup.vInfinity, setup.pericentre, kEndDistance * a));
    return setup;
}

/// A single body 0, at rest at the origin, meets a binary (bodies 1 and 2;
/// orbit a, e) whose centre of mass arrives from a population of velocity
/// dispersion sigma with a pericentre drawn uniformly in [0, maxPericentre].
/// Both orbits are oriented isotropically.  (Code 0.)
inline EncounterSetup drawSingleBinary(const Bodies& bodies, double a, double e, double maxPericentre, double sigma,
                                       Rng& rng) {
    constexpr double kSmallestPericentre = 1e-9;   // in units of a: keeps the orbit off a head-on collision
    constexpr double kPassages = 2.0;              // follow an unsettled encounter for this many round trips

    const double m1 = bodies[0].mass, m2 = bodies[1].mass, m3 = bodies[2].mass;
    EncounterSetup setup = newSetup(bodies, m2 + m3, a, e);

    const double meanAnomaly = uniform(rng, 0.0, kTwoPi);
    const Orientation innerPlane = randomOrientation(rng);
    const RelativeState inner = boundOrbit(kG * (m2 + m3), a, e, meanAnomaly, innerPlane);   // from body 1 to body 2

    setup.pericentre = std::max(uniform(rng, 0.0, maxPericentre), kSmallestPericentre * a);
    setup.vInfinity = maxwellSpeed(rng, sigma);
    const Orientation approach = randomOrientation(rng);
    const double mu = kG * (m1 + m2 + m3);
    // Never inside twice the largest pericentre, whatever radius the table calls the breakup radius.
    const double start = std::max(startDistance(a, m2 + m3, m1), 2.0 * maxPericentre);
    const Approach outer = incomingOrbit(mu, setup.vInfinity, setup.pericentre, start, approach);

    const RelativeState& centre = outer.state;   // of the binary, relative to body 0
    const double share2 = m2 / (m2 + m3), share3 = m3 / (m2 + m3);
    setBody(setup.initialState, 0, {}, {});
    setBody(setup.initialState, 1, centre.r - share3 * inner.r, centre.v - share3 * inner.v);
    setBody(setup.initialState, 2, centre.r + share2 * inner.r, centre.v + share2 * inner.v);

    setup.inclination = approach.inclination;
    // A round trip goes out to kEndDistance times the larger of a and the largest pericentre, and back.
    const double halfTrip =
        timeFromPericentre(mu, setup.vInfinity, setup.pericentre, kEndDistance * std::max(a, maxPericentre));
    setTimes(setup, mu, start, outer.timeToPericentre, (2.0 * kPassages - 1.0) * halfTrip);
    return setup;
}

// ----------------------------------------------------------------------------
// Outcome classification
// ----------------------------------------------------------------------------

/// The three ways of splitting the bodies into a pair {i, j} and a third
/// body k, in the notebooks' order: m1/m2, m1/m3, m2/m3.
struct Split {
    int i, j, k;
};
inline constexpr Split kSplits[3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};

/// Orbit of the third body of a split about the centre of mass of the pair.
struct OuterOrbit {
    Vec3 r, v;
    double mu = 0.0;       // gravitational parameter of all three masses
    double energy = 0.0;   // specific orbital energy; negative if bound
};

inline OuterOrbit outerOrbit(const Masses& m, const Vec3* r, const Vec3* v, const Split& s) {
    const double mPair = m[s.i] + m[s.j];
    OuterOrbit outer;
    outer.r = r[s.k] - (m[s.i] * r[s.i] + m[s.j] * r[s.j]) / mPair;
    outer.v = v[s.k] - (m[s.i] * v[s.i] + m[s.j] * v[s.j]) / mPair;
    outer.mu = kG * (mPair + m[s.k]);
    outer.energy = dot(outer.v, outer.v) / 2 - outer.mu / norm(outer.r);
    return outer;
}

enum class Outcome {
    BinaryM1M2,      // "m1/m2-m3": the pair m1/m2 is left bound, m3 escapes
    BinaryM1M3,      // "m1/m3-m2"
    BinaryM2M3,      // "m2/m3-m1"
    TripleM1M2,      // "(m1/m2)/m3": bound triple with inner binary m1/m2
    TripleM1M3,      // "(m1/m3)/m2"
    TripleM2M3,      // "(m2/m3)/m1"
    Ionised,         // "m1-m2-m3": all three unbound
    Unresolved       // interaction not finished at the end of the integration
};

inline const char* outcomeName(Outcome outcome) {
    static const char* names[] = {"m1/m2-m3",   "m1/m3-m2",   "m2/m3-m1", "(m1/m2)/m3",
                                  "(m1/m3)/m2", "(m2/m3)/m1", "m1-m2-m3", "UNRESOLVED / RESONANT"};
    return names[static_cast<int>(outcome)];
}

struct FinalConfiguration {
    Outcome outcome = Outcome::Unresolved;
    int pair = -1;                         // tightest bound pair, an index into kSplits; -1 if none
    double aBinary = 0.0, eBinary = 0.0;   // its orbit
    double aOuter = 0.0, eOuter = 0.0;     // third body's orbit (bound triples only)
    double mutualInclination = 0.0;        // [deg] (bound triples only)
    double recoilSpeed = 0.0;              // of the pair in the three-body centre-of-mass frame (binary + escaper only)
    double escaperDistance = 0.0;          // of the third body from the pair's centre of mass (binary + escaper only)

    bool resolved() const { return outcome != Outcome::Unresolved; }
    bool isTriple() const { return outcome >= Outcome::TripleM1M2 && outcome <= Outcome::TripleM2M3; }

    /// True if the tightest bound pair is the given split, in a binary or a triple.
    bool boundPairIs(int split) const { return resolved() && pair == split; }
};

/// Classify the end state of an encounter.  aInitial is the semi-major axis
/// of the binary before the encounter (sets the "far apart" scale when no
/// pair is bound).
inline FinalConfiguration classify(const Masses& m, const State& y, double aInitial) {
    Vec3 r[3], v[3];
    for (int i = 0; i < 3; ++i) r[i] = position(y, i), v[i] = velocity(y, i);
    const Vec3 vCom = (m[0] * v[0] + m[1] * v[1] + m[2] * v[2]) / (m[0] + m[1] + m[2]);

    PairOrbit orbit[3];
    int tightest = -1;   // bound pair with the most negative specific energy
    for (int p = 0; p < 3; ++p) {
        const Split& s = kSplits[p];
        orbit[p] = pairOrbit(r[s.i], v[s.i], m[s.i], r[s.j], v[s.j], m[s.j]);
        if (orbit[p].energy < 0 && (tightest < 0 || orbit[p].energy < orbit[tightest].energy)) tightest = p;
    }

    FinalConfiguration result;

    if (tightest < 0) {
        // No bound pair: ionised if all pairs are receding, far apart, and the
        // total energy in the centre-of-mass frame is positive.
        bool allReceding = true;
        double minSeparation = std::numeric_limits<double>::infinity();
        for (const PairOrbit& o : orbit) {
            allReceding = allReceding && dot(o.separation, o.relVelocity) > 0;
            minSeparation = std::min(minSeparation, norm(o.separation));
        }
        double energy = 0.0;
        for (int i = 0; i < 3; ++i) energy += 0.5 * m[i] * dot(v[i] - vCom, v[i] - vCom);
        for (int p = 0; p < 3; ++p) energy -= kG * m[kSplits[p].i] * m[kSplits[p].j] / norm(orbit[p].separation);
        if (energy > 0 && allReceding && minSeparation > kEscapeDistance * aInitial) result.outcome = Outcome::Ionised;
        return result;
    }

    // At least one bound pair: take the tightest as "the binary" and look at
    // the third body's orbit about its centre of mass.
    const Split& s = kSplits[tightest];
    result.pair = tightest;
    result.aBinary = orbit[tightest].semiMajorAxis;
    result.eBinary = orbit[tightest].eccentricity;

    const OuterOrbit outer = outerOrbit(m, r, v, s);
    if (outer.energy < 0) {
        const Vec3 hOut = cross(outer.r, outer.v);
        const Vec3 hIn = orbit[tightest].angMomentum;
        result.eOuter = norm(cross(outer.v, hOut) / outer.mu - outer.r / norm(outer.r));
        result.aOuter = -outer.mu / (2 * outer.energy);
        const double cosI = std::clamp(dot(hIn, hOut) / (norm(hIn) * norm(hOut)), -1.0, 1.0);
        result.mutualInclination = std::acos(cosI) * 180.0 / kPi;
        result.outcome = static_cast<Outcome>(static_cast<int>(Outcome::TripleM1M2) + tightest);
    } else if (outer.energy > 0 && dot(outer.r, outer.v) > 0 && norm(outer.r) > kEscapeDistance * result.aBinary) {
        const Vec3 vPair = (m[s.i] * v[s.i] + m[s.j] * v[s.j]) / (m[s.i] + m[s.j]);
        result.recoilSpeed = norm(vPair - vCom);
        result.escaperDistance = norm(outer.r);
        result.outcome = static_cast<Outcome>(static_cast<int>(Outcome::BinaryM1M2) + tightest);
    }
    return result;
}

/// True if nothing more can happen, so that the integration may end before
/// its time is up: three single bodies, or a binary and a body that leaves
/// it and is as far from it as an incoming body starts (startDistance, for
/// the wider of the binary left and the binary before the encounter).  From
/// there on the binary is alone, as it was before the integration began.
/// A bound triple is never settled: it is followed to the end.
inline bool settled(const EncounterSetup& setup, const FinalConfiguration& final) {
    if (final.outcome == Outcome::Ionised) return true;
    if (final.outcome > Outcome::BinaryM2M3) return false;   // a triple, or unresolved
    const Split& s = kSplits[final.pair];
    const Masses& m = setup.masses;
    const double a = std::max(final.aBinary, setup.binary.semiMajorAxis);
    return final.escaperDistance > startDistance(a, m[s.i] + m[s.j], m[s.k]);
}

/// True if the encounter left a star bound to a black hole as the tightest
/// pair, on an orbit whose pericentre lies within `factor` tidal radii: the
/// star does not survive its next passage.
inline bool starDisrupted(const Bodies& bodies, const FinalConfiguration& final, double factor) {
    if (!final.resolved() || final.pair < 0) return false;
    const Body& first = bodies[kSplits[final.pair].i];
    const Body& second = bodies[kSplits[final.pair].j];
    if (first.isStar == second.isStar) return false;   // two stars or two black holes
    return tidallyDisrupted(first.isStar ? second.mass : first.mass, final.aBinary, final.eBinary, factor);
}

// ----------------------------------------------------------------------------
// Running an encounter
// ----------------------------------------------------------------------------

struct EncounterResult {
    bool ok = true;              // the integrator reached the end time
    long steps = 0;
    int retries = 0;             // restarts with a looser tolerance (see IntegrationStats)
    double time = 0.0;           // integrated: up to setup.duration
    double initialEnergy = 0.0, finalEnergy = 0.0;
    State finalState{};
    FinalConfiguration final;    // set only if ok

    /// |dE / E| over the integration.
    double energyError() const { return ok ? std::fabs((finalEnergy - initialEnergy) / initialEnergy) : 0.0; }
};

/// Integrate an encounter and classify the result.
///
/// The outcome is looked at for the first time at setup.firstLook, and then
/// each time the time since pericentre has doubled.  The integration ends at
/// the first look that finds the encounter settled, and at setup.duration at
/// the latest.  The looks fall at times fixed by the initial conditions, so
/// every integrator stops at the same time.
///
/// With a sampler the trajectory is recorded instead, in one piece up to the
/// sampler's end time, and classified there.
template <class Integrator>
EncounterResult runEncounter(const Integrator& integrator, const EncounterSetup& setup, Sampler* sampler = nullptr) {
    constexpr double kLastStretch = 0.75;   // a look later than this share of the way to the end is the end

    EncounterResult result;
    result.finalState = setup.initialState;
    result.initialEnergy = totalEnergy(setup.masses, setup.initialState);
    if (sampler) sampler->record(setup.initialState);

    const double end = sampler ? sampler->endTime() : setup.duration;
    double sincePericentre = setup.firstLook - setup.approachTime;
    while (result.ok && result.time < end) {
        double look = setup.approachTime + sincePericentre;
        if (sampler || look - setup.approachTime > kLastStretch * (end - setup.approachTime)) look = end;

        const IntegrationStats stats =
            integrator.integrate(setup.masses, result.finalState, look - result.time, setup.binary, sampler);
        result.steps += stats.steps;
        result.retries += stats.retries;
        result.ok = stats.ok && allFinite(result.finalState);
        if (!result.ok) break;

        result.time = look;
        result.final = classify(setup.masses, result.finalState, setup.binary.semiMajorAxis);
        if (settled(setup, result.final)) break;
        sincePericentre *= 2.0;
    }
    if (result.ok) result.finalEnergy = totalEnergy(setup.masses, result.finalState);
    return result;
}

/// "Integrate this encounter": how the stages of the pipeline see the
/// integrator, whichever one it is.
using Scatter = std::function<EncounterResult(const EncounterSetup&)>;

}  // namespace bhs
