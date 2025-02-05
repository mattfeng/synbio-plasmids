def revcomp(s):
    table = str.maketrans("ATCGatcg", "TAGCtagc")

    return s.translate(table)[::-1]

modified_frag_fwd_top = "cttctagagcgtctct{outer}{inner}tgagaccggagttgac"
modified_frag_rev_top = "gttgcacgctggtctct{inner}{outer}tgagacgtactagtagcg"

unmodified_frag_fwd_top = "tgagacgtactagtagcg"
unmodified_frag_rev_top = "cttctagagcgtctct"

overhangs = {
    1: "GGAG",
    2: "TACT",
    3: "AATG",
    4: "GCTT",
    5: "CGCT"
}

outer = "AATG"
inner = "GGAG"

from_overhangs = [
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5)
]

to_overhangs = [
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5)
]

unmodified_frag_fwd = unmodified_frag_fwd_top
unmodified_frag_rev = revcomp(unmodified_frag_rev_top)
print("const_frag_fwd", unmodified_frag_fwd)
print("const_frag_rev", unmodified_frag_rev)
print()

for f5, f3 in from_overhangs:
    for t5, t3 in to_overhangs:

        outer5 = overhangs[t5]
        inner5 = overhangs[f5]
        outer3 = overhangs[t3]
        inner3 = overhangs[f3]

        if f5 == t5:
            inner5 = ""

        if f3 == t3:
            inner3 = ""

        modified_frag_fwd = modified_frag_fwd_top.format(outer=outer5, inner=inner5)
        modified_frag_rev = revcomp(modified_frag_rev_top.format(outer=outer3, inner=inner3))

        print(f"{f5}{f3}-{t5}{t3}_mod_frag_fwd", modified_frag_fwd)
        print(f"{f5}{f3}-{t5}{t3}_mod_frag_rev", modified_frag_rev)
        print()
