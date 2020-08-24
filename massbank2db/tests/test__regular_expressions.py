import re


rt_regex = '^AC\$CHROMATOGRAPHY:\s+RETENTION.*TIME\s+(\d*[.,]?\d*)\s?($|min|sec)'

rt_strs = [
    "AC$CHROMATOGRAPHY: RETENTION_TIME 23.46 min (in paper: 23.3 min)",  # Chubu_Univ/UT001973.txt
    "AC$CHROMATOGRAPHY: RETENTION_TIME 2.75",                            # AAFC/AC000338.txt
    "AC$CHROMATOGRAPHY: RETENTION_TIME 9.367 min",                       # Athens_Univ/AU592019.txt
    "AC$CHROMATOGRAPHY: RETENTION_TIME 478.2 sec",                       # BS/BS001044.txt:
]

for rt_str in rt_strs:
    print(rt_str)

    m = re.search(rt_regex, rt_str, re.IGNORECASE)
    print("\tMatch:", m)

    for i in [1, 2]:
        print("\t", m.group(i))

    print(m.groups())
