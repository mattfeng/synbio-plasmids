def revcomp(s):
    table = str.maketrans("ATCGatcg", "TAGCtagc")

    return s.translate(table)[::-1]

modified_frag_fwd_top = "cttctagagcgtctct{outer}{inner}tgagaccggagttgac"
modified_frag_rev_top = "gttgcacgctggtctct{inner}{outer}tgagacgtactagtagcg"

unmodified_frag_fwd_top = "tgagacgtactagtagcg"
unmodified_frag_rev_top = "cttctagagcgtctct"

outer = "AATG"
inner = "GGAG"

modified_frag_fwd = modified_frag_fwd_top.format(outer=outer, inner=inner)
modified_frag_rev = revcomp(modified_frag_rev_top.format(outer=outer, inner=inner))

unmodified_frag_fwd = unmodified_frag_fwd_top
unmodified_frag_rev = revcomp(unmodified_frag_rev_top)

print("modified_frag_fwd", modified_frag_fwd)
print("modified_frag_rev", modified_frag_rev)
print("unmodified_frag_fwd", unmodified_frag_fwd)
print("unmodified_frag_rev", unmodified_frag_rev)
