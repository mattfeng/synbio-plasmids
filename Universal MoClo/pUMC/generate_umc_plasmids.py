# unfinished

from Bio import SeqIO
from Bio.Seq import Seq

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


snapgene_file = "4bp/4bp_A.dna"

record = SeqIO.read(snapgene_file, "snapgene")

print(f"Record ID: {record.id}")
print(f"Sequence: {record.seq}")
print(f"Description: {record.description}")
print(f"Length: {len(record.seq)}")

outer_5_idx = 76
outer_3_idx = 923
inner_5_idx = 80
inner_3_idx = 919

overhangs = {
    1: "GGAG",
    2: "TACT",
    3: "AATG",
    4: "GCTT",
    5: "CGCT",
}

def all_pairs(n):
    for i in range(1, n):
        for j in range(i + 1, n):
            yield (i, j)

def delete_range_from_seqrecord(seq_record: SeqRecord, start: int, end: int) -> SeqRecord:
    """
    Deletes a contiguous range of bases from a SeqRecord, adjusts annotations:
    - Removes features entirely within the range.
    - Truncates features partially overlapping the range.
    - Shifts features after the range.

    Args:
        seq_record (SeqRecord): The input sequence record.
        start (int): Start of the range to delete (0-based, inclusive).
        end (int): End of the range to delete (0-based, exclusive).

    Returns:
        SeqRecord: A new SeqRecord with the range deleted and annotations updated.
    """
    # Create a new sequence excluding the specified range
    new_seq = seq_record.seq[:start] + seq_record.seq[end:]

    # Adjust features
    new_features = []
    for feature in seq_record.features:
        feature_start = feature.location.start
        feature_end = feature.location.end

        # If the feature is completely outside the deletion range, adjust its position
        if feature_end <= start:
            # Feature before the deleted region remains unchanged
            new_features.append(feature)
        elif feature_start >= end:
            # Feature after the deleted region shifts by the deletion length
            shift = end - start
            new_location = FeatureLocation(
                feature_start - shift,
                feature_end - shift,
                strand=feature.location.strand
            )
            new_features.append(SeqFeature(location=new_location, type=feature.type, qualifiers=feature.qualifiers))
        elif feature_start < start and feature_end > end:
            # Feature overlaps the range, truncate the deleted portion
            new_location = FeatureLocation(
                feature_start,
                feature_end - (end - start),
                strand=feature.location.strand
            )
            new_features.append(SeqFeature(location=new_location, type=feature.type, qualifiers=feature.qualifiers))
        # Features completely within the range are not added

    # Create a new SeqRecord with the updated sequence and features
    new_record = SeqRecord(
        seq=new_seq,
        id=seq_record.id,
        name=seq_record.name,
        description=seq_record.description,
        annotations=seq_record.annotations
    )
    new_record.features = new_features

    return new_record

for out_type in all_pairs(6):
    for in_type in all_pairs(6):

        out_5_type, out_3_type = out_type
        in_5_type, in_3_type = in_type

        outer_5 = overhangs[out_5_type]
        outer_3 = overhangs[out_3_type]
        inner_5 = overhangs[in_5_type]
        inner_3 = overhangs[in_3_type]

        seq = list(record.seq)

        seq[outer_5_idx - 1:outer_5_idx + 3] = list(outer_5)
        seq[outer_3_idx - 1:outer_3_idx + 3] = list(outer_3)
        seq[inner_5_idx - 1:inner_5_idx + 3] = list(inner_5)
        seq[inner_3_idx - 1:inner_3_idx + 3] = list(inner_3)

        record.seq = Seq("".join(seq))

        # remove duplicate scars if possible
        if out_5_type == in_5_type:
            record = delete_range_from_seqrecord(
                record,
                inner_5_idx - 1,
                inner_5_idx + 3
                )
        if out_3_type == in_3_type:
            record = delete_range_from_seqrecord(
                record,
                inner_3_idx - 1,
                inner_3_idx + 3
                )

        SeqIO.write(record, f"4bp/A/pUMC-A{in_5_type}{in_3_type}-B{out_5_type}{out_3_type}.gbk", "genbank")
