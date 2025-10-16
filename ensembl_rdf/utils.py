import urllib.parse
import sys
from datetime import datetime

def quote_str(string):
    # Enclose a string with double-quotation marks.
    # Double-quotation marks within the input string are escaped with backslashes.
    return "\"" + string.replace("\"", "\\\"") + "\""


def percent_encode(string):
    #return re.sub(r"([()])", r"\\\1", string)
    return urllib.parse.quote(string, safe="")


def strand2faldo(s):
    if s == "1":
        return "faldo:ForwardStrandPosition"
    elif s == "-1":
        return "faldo:ReverseStrandPosition"
    else:
        print(f"Error: Invalid argument \"{s}\" for strand2faldo", file=sys.stderr)
        sys.exit(1)


def log_time(text):
    dt_now = datetime.now()
    print(f"[{dt_now}] {text}", file=sys.stderr)
    return


class Bnode:
    def __init__(self):
        self.properties = []

    def add(self, tpl):
        # `tpl` is a tuple of strings (e.g. ("rdf:type", "owl:Class"))
        self.properties.append(tpl)

    def serialize(self, level=1):
        s = "[\n"
        indent = "    " * level
        for tpl in self.properties:
            s += indent + tpl[0] + " " + tpl[1] + " ;\n"
        s += "    " * (level-1) + "]"
        return s
