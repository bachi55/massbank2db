def get_precursor_mz(exact_mass, precursor_type):
    """
    Calculate precursor mz based on exact mass and precursor type

    :param exact_mass: float, exact mass of compound of interest

    :param precursor_type: string, precursor type, e.g. [M+H]+

    :return: scalar, neutral mass of compound
    """
    # Values are taken from:
    # https://docs.google.com/spreadsheets/d/1r4dPw1shIEy_W2BkfgPsihinwg-Nah654VlNTn8Gxo0/edit?usp=sharing
    d = {'[M+H]+': 1.007276,
         '[M-H]-': -1.00728,
         '[M+HCOO]-': -44.9982,
         '[M+CH3COO]-': -59.01385,
         '[M+Na]+': 22.98922,
         '[M+CH3COOH-H]-': None,
         '[M]+': 0,
         '[M-2H]-': -1.00728,
         '[M-H2O+H]+': None,
         '[M+NH4]+': 18.03383,
         '[M-2H2O+H]+': None,
         '[M-H-C6H10O5]-': None,
         '[M-H+CH2O2]-': None,
         '[M-H+C2H2O]-': None,
         '[M-CH3]-': None,
         '[M+H-H2O]+': None,
         '[M+H-NH3]+': None,
         '[M+H-C9H10O5]+': None,
         '[M-C6H10O5+H]+': None
    }

    try:
        return exact_mass + d[precursor_type]
    except KeyError as e:
        print(e)
        return False
