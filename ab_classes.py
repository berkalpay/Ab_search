from antibodies import Gene, Class

classes = [
    Class(
        name="VH1-18 QXXV",
        vs=[Gene("IGHV1-18")],
        cdr_length=(18, 22),
        cdr_signature=(98, "Q..(V|I)(.*)"),
    ),
]
