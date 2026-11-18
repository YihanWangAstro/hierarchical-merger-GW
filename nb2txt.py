#!/usr/bin/env python3
"""Turn a Mathematica notebook (.nb) into text that can be read without Mathematica.

    python nb2txt.py "2 Code 0.nb" "2 Code 0.txt"
    python nb2txt.py notebook.nb                      # to the terminal
    python nb2txt.py notebook.nb out.txt --max-output 20000

Input cells are written in full, as the code one would type. Output and Print
cells are cut after --max-output characters (default 4000); graphics are only
named.
"""
import re
import sys
import unicodedata

NAMED = {
    "IndentingNewLine": "\n", "InvisibleSpace": "", "InvisibleTimes": " ", "InvisibleApplication": "@",
    "InvisibleComma": ",", "NoBreak": "", "Null": "", "LineSeparator": "\n",
    "Rule": "->", "RuleDelayed": ":>", "LessEqual": "<=", "GreaterEqual": ">=", "NotEqual": "!=",
    "Equal": "==", "LongEqual": "==", "And": "&&", "Or": "||", "Not": "!", "Times": "*", "Cross": "×",
    "Function": "|->", "LeftDoubleBracket": "[[", "RightDoubleBracket": "]]", "Infinity": "Infinity",
    "Degree": "Degree", "ExponentialE": "E", "ImaginaryI": "I", "ImaginaryJ": "I", "Pi": "π",
    "DifferentialD": "d", "Integral": "∫", "Sum": "∑", "Product": "∏", "PartialD": "∂", "Sqrt": "√",
    "Element": "∈", "Bullet": "•", "Ellipsis": "…", "CenterDot": "·", "Prime": "′", "Micro": "µ",
    "Alpha": "α", "Beta": "β", "Gamma": "γ", "Delta": "δ", "Epsilon": "ϵ", "CurlyEpsilon": "ε", "Zeta": "ζ",
    "Eta": "η", "Theta": "θ", "CurlyTheta": "ϑ", "Iota": "ι", "Kappa": "κ", "Lambda": "λ", "Mu": "μ",
    "Nu": "ν", "Xi": "ξ", "Omicron": "ο", "Rho": "ρ", "Sigma": "σ", "FinalSigma": "ς", "Tau": "τ",
    "Upsilon": "υ", "Phi": "ϕ", "CurlyPhi": "φ", "Chi": "χ", "Psi": "ψ", "Omega": "ω",
    "CapitalAlpha": "Α", "CapitalBeta": "Β", "CapitalGamma": "Γ", "CapitalDelta": "Δ", "CapitalTheta": "Θ",
    "CapitalLambda": "Λ", "CapitalXi": "Ξ", "CapitalPi": "Π", "CapitalSigma": "Σ", "CapitalPhi": "Φ",
    "CapitalPsi": "Ψ", "CapitalOmega": "Ω", "HBar": "ℏ", "Angstrom": "Å", "PlusMinus": "±", "Proportional": "∝",
    "TildeTilde": "≈", "Tilde": "~", "RightArrow": "→", "LeftArrow": "←", "Minus": "-", "LessLess": "<<",
    "GreaterGreater": ">>", "Conjugate": "*", "Transpose": "^T", "LeftBracketingBar": "|",
    "RightBracketingBar": "|", "VerticalSeparator": "|", "Dash": "-", "LongDash": "—", "Odot": "⊙", "CirclePlus": "⊕",
    "CircleDot": "⊙", "Sun": "☉", "RawEscape": "", "AliasDelimiter": "", "SpanFromLeft": "", "Placeholder": "□",
    "LeftAssociation": "<|", "RightAssociation": "|>", "DoubleStruckCapitalR": "ℝ", "ScriptL": "ℓ",
    "ScriptCapitalM": "ℳ", "ScriptCapitalR": "ℛ", "ScriptCapitalE": "ℰ", "DotlessI": "ı", "Star": "⋆",
    "Earth": "⊕", "Square": "□", "Del": "∇", "Distributed": "~", "Conditioned": "|", "Divide": "/",
}

TOKEN = re.compile(
    r"""\s*(?:
      (?P<comment>\(\*.*?\*\))
    | (?P<string>"(?:[^"\\]|\\.)*")
    | (?P<number>[+-]?(?:\d+\.?\d*|\.\d+)(?:`{1,2}[+-]?\d*\.?\d*)?(?:\*\^[+-]?\d+)?)
    | (?P<symbol>[A-Za-z$][A-Za-z0-9$`]*)
    | (?P<op>->|:>|\[|\]|\{|\}|,|.)
    )""",
    re.S | re.X,
)


class Parser:
    def __init__(self, text):
        self.toks = []
        for m in TOKEN.finditer(text):
            kind = m.lastgroup
            if kind and kind != "comment":
                self.toks.append((kind, m.group(kind)))
        self.i = 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else ("eof", "")

    def take(self):
        tok = self.peek()
        self.i += 1
        return tok

    def expr(self):
        left = self.term()
        while self.peek()[1] in ("->", ":>"):
            op = self.take()[1]
            right = self.expr()
            left = ("Rule" if op == "->" else "RuleDelayed", [left, right])
        return left

    def args(self, closer):
        out = []
        if self.peek()[1] == closer:
            self.take()
            return out
        while True:
            out.append(self.expr())
            tok = self.take()
            if tok[1] == closer:
                return out
            if tok[1] != ",":
                # tolerate arithmetic inside option values: swallow until , or closer
                depth = 0
                while True:
                    if tok[0] == "eof":
                        raise SyntaxError("unexpected end")
                    if tok[1] in "[{(":
                        depth += 1
                    elif tok[1] in "]})":
                        if depth == 0 and tok[1] == closer:
                            return out
                        depth -= 1
                    elif tok[1] == "," and depth == 0:
                        break
                    tok = self.take()

    def term(self):
        kind, val = self.take()
        if kind == "string":
            node = ("str", val[1:-1])
        elif kind == "number":
            node = ("num", val)
        elif kind == "symbol":
            node = ("sym", val)
        elif val == "{":
            node = ("List", self.args("}"))
        else:
            node = ("sym", val)
        while self.peek()[1] == "[":
            self.take()
            head = node[1] if node[0] == "sym" else node
            node = (head, self.args("]"))
        return node


def unescape(s):
    s = s.replace("\\\n", "")
    out, i = [], 0
    while i < len(s):
        ch = s[i]
        if ch != "\\":
            out.append(ch)
            i += 1
            continue
        nxt = s[i + 1] if i + 1 < len(s) else ""
        if nxt == "[":
            j = s.index("]", i)
            name = s[i + 2:j]
            if name in NAMED:
                out.append(NAMED[name])
            else:
                out.append("\\[" + name + "]")
            i = j + 1
        elif nxt == ":":
            out.append(chr(int(s[i + 2:i + 6], 16)))
            i += 6
        elif nxt in "<>":
            i += 2
        elif nxt == "n":
            out.append("\n")
            i += 2
        elif nxt == "t":
            out.append("\t")
            i += 2
        elif nxt in "!*()":      # linear-syntax markers \!\(\*...\)
            i += 2
        else:
            out.append(nxt)
            i += 2
    return "".join(out)


SIMPLE = re.compile(r"^[\w.`$ϵ]+$|^-?[\d.]+$")


def paren(s):
    return s if SIMPLE.match(s) else "(" + s + ")"


def number(s):
    s = re.sub(r"`{1,2}[+-]?\d*\.?\d*", "", s)
    return s.replace("*^", "e")


def render(node):
    if node is None:
        return ""
    kind, val = node
    if kind == "str":
        text = unescape(val)
        if re.fullmatch(r"[+-]?(\d+\.?\d*|\.\d+)`{1,2}[-+\d.]*(\*\^[+-]?\d+)?", text):
            return number(text)
        return text
    if kind == "num":
        return number(val)
    if kind == "sym":
        return val
    args = val
    r = lambda k: render(args[k]) if k < len(args) else ""
    if kind == "List":
        return "{" + ",".join(render(a) for a in args) + "}"
    if kind == "RowBox":
        inner = args[0]
        if inner[0] == "List":
            return "".join(render(a) for a in inner[1])
        return render(inner)
    if kind == "BoxData":
        inner = args[0]
        if inner[0] == "List":
            return "".join(render(a) for a in inner[1])
        return render(inner)
    if kind == "SuperscriptBox":
        base, exp = r(0), r(1)
        if exp == "′":
            return base + "'"
        return paren(base) + "^" + paren(exp)
    if kind == "SubscriptBox":
        return "Subscript[" + r(0) + "," + r(1) + "]"
    if kind == "SubsuperscriptBox":
        return "Subsuperscript[" + r(0) + "," + r(1) + "," + r(2) + "]"
    if kind == "FractionBox":
        return "((" + r(0) + ")/(" + r(1) + "))"
    if kind == "SqrtBox":
        return "Sqrt[" + r(0) + "]"
    if kind == "RadicalBox":
        return "(" + r(0) + ")^(1/(" + r(1) + "))"
    if kind in ("UnderoverscriptBox", "UnderscriptBox", "OverscriptBox"):
        return kind[:-3] + "[" + ",".join(render(a) for a in args[:3] if a[0] not in ("Rule", "RuleDelayed")) + "]"
    if kind in ("StyleBox", "TagBox", "InterpretationBox", "FormBox", "PaneBox", "ItemBox", "FrameBox",
                "AdjustmentBox", "TooltipBox", "ErrorBox", "PanelBox"):
        return r(0)
    if kind == "GridBox":
        rows = args[0][1] if args and args[0][0] == "List" else []
        return "{" + ",".join("{" + ",".join(render(c) for c in row[1]) + "}" for row in rows if row[0] == "List") + "}"
    if kind == "TemplateBox":
        tag = render(args[1]) if len(args) > 1 else ""
        items = args[0][1] if args and args[0][0] == "List" else []
        return "Template<" + tag + ">[" + ",".join(render(a) for a in items) + "]"
    if kind in ("Rule", "RuleDelayed"):
        return r(0) + ("->" if kind == "Rule" else ":>") + r(1)
    if isinstance(kind, str) and kind.endswith("Box") or kind in ("DynamicModuleBox", "GraphicsBox", "Graphics3DBox"):
        return "<" + kind + ">"
    head = kind if isinstance(kind, str) else render(kind)
    return head + "[" + ",".join(render(a) for a in args) + "]"


def match_bracket(text, start):
    """Index just past the ] matching the [ at text[start]."""
    depth, i, n = 0, start, len(text)
    while i < n:
        ch = text[i]
        if ch == '"':
            i += 1
            while text[i] != '"':
                i += 2 if text[i] == "\\" else 1
        elif ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise SyntaxError("unbalanced")


def cells(text):
    pos = 0
    pattern = re.compile(r"Cell\[(?!CellGroupData)")
    while True:
        m = pattern.search(text, pos)
        if not m:
            return
        end = match_bracket(text, m.end() - 1)
        yield text[m.start():end]
        pos = end


def main():
    argv = sys.argv[1:]
    max_output = 4000   # characters of an output cell that are kept
    if "--max-output" in argv:
        at = argv.index("--max-output")
        max_output = int(argv[at + 1])
        del argv[at:at + 2]
    if not argv:
        sys.exit(__doc__)
    text = open(argv[0], encoding="utf-8", errors="replace").read()
    out = open(argv[1], "w", encoding="utf-8") if len(argv) > 1 else sys.stdout

    for cell in cells(text):
        body_end = match_bracket(cell, cell.index("[", 5)) if cell[5:12] in ("BoxData", "TextDat") else None
        if body_end is None:
            m = re.match(r'Cell\[("(?:[^"\\]|\\.)*")\s*,\s*"(\w+)"', cell, re.S)
            if m:
                out.write("===== %s =====\n%s\n\n" % (m.group(2), unescape(m.group(1)[1:-1])))
            continue
        rest = cell[body_end:]
        sm = re.match(r'\s*,\s*"(\w+)"', rest)
        style = sm.group(1) if sm else "?"
        lm = re.search(r'CellLabel->\s*"([^"]*)"', rest)
        label = lm.group(1) if lm else ""
        body = cell[5:body_end]
        header = "===== %s %s =====\n" % (style, label)
        if style not in ("Input", "Code") and len(body) > max_output * 12:
            kinds = sorted(set(re.findall(r"\b(\w*Graphics\w*Box|GridBox|DynamicModuleBox)\b", body)))
            out.write(header + "<%d chars of boxes: %s>\n\n" % (len(body), ", ".join(kinds) or "data"))
            continue
        try:
            rendered = render(Parser(body).expr())
        except Exception as err:  # noqa: BLE001
            out.write(header + "<could not render: %s; %d chars>\n\n" % (err, len(body)))
            continue
        rendered = unicodedata.normalize("NFC", rendered)
        if style not in ("Input", "Code") and len(rendered) > max_output:
            rendered = rendered[:max_output] + " ...<%d chars total>" % len(rendered)
        out.write(header + rendered.strip("\n") + "\n\n")


if __name__ == "__main__":
    main()
