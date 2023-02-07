waters=(
    {
        "O":14,
        "H1":15,
        "H2":16,
    }, # Away from cage
    {
        "O":17,
        "H1":18,
        "H2":19,
    }, # Hydrogen bonded to amine
    {
        "O":143,
        "H1":144,
        "H2":145,
    }, # Hydrogen bonded to phenol
    {
        "O":146,
        "H1":147,
        "H2":148,
    }, # Showed significant diffusion
    {
        "O":272,
        "H1":273,
        "H2":274,
    }, # Hydrogen bonded to phenol
    {
        "O":275,
        "H1":276,
        "H2":277,
    }, # Showed significant diffusion
    {
        "O":401,
        "H1":402,
        "H2":403,
    }, # Away from cage
    {
        "O":404,
        "H1":405,
        "H2":406,
    }, # Showed significant diffusion
)

TFAs = (
    # Two free TFA
    {
        "F1": 0,
        "F2": 1,
        "F3": 2,
        "O1": 3,
        "O2": 4,
        "C1": 5,
        "C2": 6,
    }, # Free

    {
        "F1": 7,
        "F2": 8,
        "F3": 9,
        "O1": 10,
        "O2": 11,
        "C1": 12,
        "C2": 13,
    }, # Bound

    {
        "F1": 129,
        "F2": 130,
        "F3": 131,
        "O1": 132,
        "O2": 133,
        "C1": 134,
        "C2": 135,
    }, # Bound

    {
        "F1": 136,
        "F2": 137,
        "F3": 138,
        "O1": 139,
        "O2": 140,
        "C1": 141,
        "C2": 142,
    }, # Bound

    {
        "F1": 258,
        "F2": 259,
        "F3": 260,
        "O1": 261,
        "O2": 262,
        "C1": 263,
        "C2": 264,
    }, # Bound

    {
        "F1": 265,
        "F2": 266,
        "F3": 267,
        "O1": 268,
        "O2": 269,
        "C1": 270,
        "C2": 271,
    }, # Bound

    {
        "F1": 387,
        "F2": 388,
        "F3": 389,
        "O1": 390,
        "O2": 391,
        "C1": 392,
        "C2": 393,
    }, # Free

    {
        "F1": 394,
        "F2": 395,
        "F3": 396,
        "O1": 397,
        "O2": 398,
        "C1": 399,
        "C2": 400,
    }, # Bound
)

phenols = (
    {"O":24, "H":25}, #
    {"O":411, "H":412}, #
    {"O":153, "H":154}, # H-bonded to water O 272
    {"O":282, "H":283}, # H-bonded to water O 143  
)

amines = (
    {"N":410, "H1":284, "H2":285}, #H1 H-bond to TFA O 397
    {"N":409, "H":128}, # H-bond to TFA O 261
    {"N":407, "H1":28, "H2":29}, #H1 H-bond to TFA O 261, H2 H-bonded to TFA O 269
    {"N":152, "H1":26, "H2":27}, #
    {"N":149, "H1":286, "H2":287}, #
    {"N":151, "H":386,}, #

    {"N":280, "H":257}, #
    {"N":278, "H1":157, "H2":158}, #
    {"N":281, "H1":413, "H2":414}, #
    
    {"N":23, "H1":155, "H2":156}, # H1 155 H-bonded to TFA O 10 H2 156 H-bonded to water O 17
    {"N":22, "H1":515}, # H 149 H-bonded to TFA O 132
    {"N":20, "H1":415, "H2":416}, # H1 415 H-bonded to TFA O 132 H2 416 H-bonded to TFA O 140
)