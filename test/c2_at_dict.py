waters=(
    {
        "O":597,
        "H1":598,
        "H2":599,
    }, # Hydrogen bonded to TFA-O donor
    {
        "O":594,
        "H1":595,
        "H2":596,
    }, # Hydrogen bonded to TFA-O donor, hydrogen bonded to phenol acceptor
    {
        "O":302,
        "H1":303,
        "H2":304,
    }, # Hydrogen bonded to two TFA-Os donor, hydrogen bonded to phenol acceptor
    {
        "O":607,
        "H1":608,
        "H2":609,
    }, # Hydrogen bonded to two TFA-Os donor, hydrogen bonded to phenol acceptor
    {
        "O":289,
        "H1":290,
        "H2":291,
    }, # Hydrogen bonded to TFA-O donor, hydrogen bonded to phenol acceptor
    {
        "O":292,
        "H1":293,
        "H2":294,
    }, # Hydrogen bonded to TFA-O donor
)

TFAs = (
    # Four free TFAs 3, 6, 7, 8 corresponding to No 4 TFA
    {
        "F1": 282,
        "F2": 284,
        "F3": 285,
        "O1": 286,
        "O2": 288,
        "C1": 283,
        "C2": 287,
    }, # Three H-bonds, O2 H-bonded to water O 594, O1 H-bonded to water O 302 and amine N 256

    {
        "F1": 276,
        "F2": 278,
        "F3": 279,
        "O1": 275,
        "O2": 281,
        "C1": 277,
        "C2": 280,
    }, # Bound

    {
        "F1": 295,
        "F2": 296,
        "F3": 297,
        "O1": 298,
        "O2": 299,
        "C1": 300,
        "C2": 301,
    }, # No. 4 free

    {
        "F1": 1,
        "F2": 2,
        "F3": 3,
        "O1": 0,
        "O2": 4,
        "C1": 5,
        "C2": 6,
    }, # No. 1 bound to amine

    {
        "F1": 269,
        "F2": 272,
        "F3": 274,
        "O1": 268,
        "O2": 270,
        "C1": 271,
        "C2": 273,
    }, # Bound to amine, low rotate density

    {
        "F1": 7,
        "F2": 8,
        "F3": 9,
        "O1": 10,
        "O2": 13,
        "C1": 11,
        "C2": 12,
    }, # No. 4 Free

    {
        "F1": 600,
        "F2": 601,
        "F3": 602,
        "O1": 603,
        "O2": 604,
        "C1": 605,
        "C2": 606,
    }, # No. 4 Free

    {
        "F1": 312,
        "F2": 313,
        "F3": 314,
        "O1": 315,
        "O2": 318,
        "C1": 316,
        "C2": 317,
    }, # No. 4 Free

    {
        "F1": 574,
        "F2": 577,
        "F3": 579,
        "O1": 573,
        "O2": 575,
        "C1": 576,
        "C2": 578,
    }, # No. 1 bound to amine
    
    {
        "F1": 587,
        "F2": 589,
        "F3": 590,
        "O1": 591,
        "O2": 593,
        "C1": 588,
        "C2": 592,
    }, # bound
    
    {
        "F1": 581,
        "F2": 583,
        "F3": 584,
        "O1": 580,
        "O2": 586,
        "C1": 582,
        "C2": 585,
    },  # bound 
    
    {
        "F1": 306,
        "F2": 307,
        "F3": 308,
        "O1": 305,
        "O2": 309,
        "C1": 310,
        "C2": 311,
    }, # No. 3 bound to amine 333
)

phenols = (
    {"O":31, "H":32}, # H-bonded to water O 302
    {"O":458, "H":459}, # H-bonded to water O 594
    {"O":153, "H":154}, # H-bonded to water O 289
    {"O":336, "H":337}, # H-bonded to water O 607    
)

TFA_fracs = (
    {"F": 329,},
    {"F": 24,},
    {"F1": 321,"F2": 322,"F3": 323,},
    {"F1": 16,"F2": 17,"F3": 18,},
    {"F1": 26,"F2": 27,"O1": 21, "O2":22}, # O Hbonded to amine
    {"F1": 331,"F2": 332,"O1": 326, "O2":327}, # O Hbonded to amine
    {"O1": 320, "O2": 319}, # O1 H-bonded to TFA N 347, O2 H-bonded to N 333
    {"O1": 14, "O2": 15}, # O1 H-bonded to TFA N 28, H bonded amine N 42
)

amines = (
    {"N":347, "H":348}, #H-bond to TFA O 320
    {"N":338, "H":339}, #
    {"N":333, "H1":334, "H2":335}, #H1 H-bond to TFA O 305, H2 H-bonded to TFA O 319
    {"N":82, "H1":83, "H2":84}, #
    {"N":561, "H1":562, "H2":563}, # H 563 H-bonded to TFA O 591
    {"N":456, "H":457,}, # H 457 H-bonded to TFA O 580
    {"N":468, "H1":469, "H2":470}, # H 470 H-bonded to TFA O 327
    {"N":453, "H1":454, "H2":455}, # H 454 H-bonded to TFA O 326

    {"N":42, "H":43}, #H-bond to TFA O 15
    {"N":33, "H":34}, #
    {"N":28, "H1":29, "H2":30}, #H1 H-bond to TFA O 0, H2 H-bonded to TFA O 14
    {"N":151, "H":152}, #H-bond to TFA O 275
    {"N":387, "H1":388, "H2":389}, #
    
    {"N":163, "H1":164, "H2":165}, # H 165 H-bonded to TFA O 22
    {"N":148, "H1":149, "H2":150}, # H 149 H-bonded to TFA O 21
    {"N":256, "H1":257, "H2":258}, # H 258 H-bonded to TFA O 286
)