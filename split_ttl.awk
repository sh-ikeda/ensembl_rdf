### usage: awk -f split_ttl.awk [-v limit=NUM] INPUT_TTL

BEGIN {
    if (!limit)
        limit = 10000000
    to = 0
}

/^@prefix/ {
    prefixes = prefixes $0 "\n"
    next
}

{
    n++
}

n>to && !/^[\] ]/ {
    to += limit
    filenum++
    output_filename = gensub(/ttl$/, "", "g", FILENAME) filenum ".ttl"
    print "Output to " output_filename > "/dev/stderr"
    print prefixes > output_filename
}

{
    print > output_filename
}
