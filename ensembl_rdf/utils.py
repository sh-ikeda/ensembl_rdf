import urllib.parse
import sys

def quote_str(string):
    # Enclose a string with double-quotation marks.
    # Double-quotation marks within the input string are escaped with backslashes.
    return "\"" + string.replace("\"", "\\\"") + "\""


def percent_encode(string):
    #return re.sub(r"([()])", r"\\\1", string)
    return urllib.parse.quote(string)


def strand2faldo(s):
    if s == "1":
        return "faldo:ForwardStrandPosition"
    elif s == "-1":
        return "faldo:ReverseStrandPosition"
    else:
        print(f"Error: Invalid argument \"{s}\" for strand2faldo", file=sys.stderr)
        sys.exit(1)


