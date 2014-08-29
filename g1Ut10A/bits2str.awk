function bits2str(bits)
{
    if (bits == 0)
        return "0"
    mask = 1
    for (; bits != 0; bits = rshift(bits, 1))
        data = (and(bits, mask) ? "1" : "0") data
    while ((length(data) % 8) != 0)
        data = "0" data
    return data
}
{
    printf("%s %s\n", $1, bits2str($1))
}