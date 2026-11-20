// Units, physical constants and model parameters (values from the notebooks).
#pragma once

namespace bhs {

// ---- unit system: AU, yr, Msun ---------------------------------------------
inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kTwoPi = 2.0 * kPi;
inline constexpr double kG = 4.0 * kPi * kPi;           // gravitational constant
inline constexpr double kSpeedOfLight = 6.32397263e4;   // AU / yr
inline constexpr double kParsec = 206264.8;             // AU

// ---- nuclear star cluster --------------------------------------------------
inline constexpr double kStarMass = 1.0;                // field and companion stars
inline constexpr double kStarRadius = 0.005;            // their radius [AU]: "RS" of the notebooks, 1.075 solar radii
inline constexpr double kSmbhMass = 4.0e6;              // central black hole
inline constexpr double kCuspSlope = 1.75;              // stars:       n ~ r^-gamma
inline constexpr double kBhCuspSlope = kCuspSlope;      // black holes: same slope
inline constexpr double kBhNumberFraction = 0.01;       // black holes per star
inline constexpr double kDensityNorm = 1.35e6;          // stars per pc^3 ...
inline constexpr double kDensityNormRadius = 0.25;      // ... at this radius [pc]
inline constexpr double kMinFieldBhMass = 5.0;          // 1G black holes: dN/dm ~ m^-2
inline constexpr double kMaxFieldBhMass = 45.0;         //   between these two masses
inline constexpr double kRemnantMassFraction = 0.95;    // of two merging black holes, left in the remnant

// ---- radial grid -------------------------------------------------------------
// Ten radii [AU], log-spaced from 0.01 to 2 pc: Routs of the notebooks,
// Exp[Subdivide[Log[0.01 PC], Log[2 PC], 9]], as the numbers Mathematica gives.
// They are written out because exp and log round differently from one
// machine to the next.  Tidal capture (tc) is worked out at all ten; the
// outer six are the radii of Code 0, as its input table holds them.
inline constexpr double kGridRadii[] = {2062.648,           3716.166119662855,  6695.224114308461, 12062.438679378582,
                                        21732.271304079746, 39153.908143096894, 70541.56932921596, 127091.09356956737,
                                        228973.44386154477, 412529.59999999986};

// ---- three-body encounters (lengths in units of the binary's a) ------------
inline constexpr double kStartTidalRadii = 10.0;        // the incoming body starts this many tidal radii of the
                                                        //   binary away (see startDistance), not this many a
inline constexpr double kFirstLook = 2.0;               // the outcome is first looked at when a body leaving as
                                                        //   it came is this many start distances away
inline constexpr double kEndDistance = 100.0;           // an encounter not settled by then ends when a body
                                                        //   leaving as it came is this far away
inline constexpr double kMinPericentre = 1e-4;          // pericentre drawn uniformly in
inline constexpr double kMaxPericentre = 2.0;           //   [kMinPericentre, kMaxPericentre]
inline constexpr double kEscapeDistance = 10.0;         // a body this far away may count as escaped

}  // namespace bhs
