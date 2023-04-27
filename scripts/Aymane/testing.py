import re


for line in genbank:
    if re.match("^LOCUS", line):
        locus = re.sub("^LOCUS\s*", "", line)
        print(locus)
    elif re.match("^ACCESSION", line):
        accession = re.sub("^ACCESSION\s*", "", line)
    elif re.match("^REFERENCE", line):
        reference = re.sub("^REFERENCE\s*", "", line, flags=re.MULTILINE)
        print(reference)
    elif re.match("^FEATURES", line):
        fields = parse_annotation(annotation)
        features = parse_features(fields["FEATURES"])
        for feature in features:
            featurename = re.search(r"^{5}(\S+)", feature).group(1)
            print(feature)
            if re.search("keywords", locus, re.IGNORECASE) or re.search("keywords", reference, re.IGNORECASE) or re.search("keywords", feature, re.IGNORECASE):
                print(accession)

def parse_features(features):
    result = []
    while re.search(r"^{5}\S.*\n(^{21}\S.*\n)*", features, re.MULTILINE):
        feature = re.search(r"^{5}\S.*\n(^{21}\S.*\n)*", features, re.MULTILINE).group()
        result.append(feature)
        features = re.sub(feature, "", features, count=1)
    return result
