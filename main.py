
import re
import math

class Atom:
    def __init__(self, data):
        self.symbol:str = data.get("name")
        self.atomic_number :int = data.get("atomic_number", 0)
        self.electrons : int = data.get("electrons", self.atomic_number)
        self.charge: int = data.get("charge", self.atomic_number - self.electrons)
        self.c_s: str = data.get("a")
        self.a_s: list[str] = data.get("c")
        self.b_s: list[str] = data.get("d")

        self.element: str = data.get("symbol", self.symbol)
        self.c: list[float] = []
        self.a: list[float] = []
        self.b: list[float] = []
        self.source: str = ""
        self.upper_limit: float = -1
        self.lower_limit: float = 0

        self._clean_name_symbol()
        self._make_numeric()

    def _make_numeric(self):
        self.c = float(self.c_s)
        self.a = [float(c) for c in self.a_s]
        self.b = [float(d) for d in self.b_s]

    def _clean_name_symbol(self):
        first_bracket = self.symbol.find("(")
        self.element = self.symbol[0:first_bracket]
        self.symbol = self.symbol.replace("(","")
        self.symbol = self.symbol.replace(")","")

    def f_s(self, s)-> float:
        r_B = 0.529177210544
        gaussians = self.c + sum([a*math.exp((-b * s**2)) for a,b in zip(self.a, self.b)])
        return self.electrons - 8 * math.pi * r_B * s**2 * gaussians

    def f_d(self, d) -> float:
        return self.f_s(1/(2 * d))

    def __str__(self):
        return f"{self.symbol}|{self.atomic_number}|{self.electrons}|{self.charge}|{self.element}|{self.c_s}|{len(self.a_s)}:{self.a_s}|{len(self.b_s)}:{self.b_s}|{self.source}"

    def cif_line(self) -> str:
        return f"{self.symbol:<4}  {self.c_s}  [ {'  '.join(self.a_s)} ]  [ {'  '.join(self.b_s)} ]  {self.lower_limit}  {self.upper_limit}  '{self.source}'"


def atom_stream_neutral(lines):
    current = {}
    current_key = None
    first_key_of_record = True

    for raw_line in lines:
        line = raw_line.strip()

        if not line:
            if current:
                yield current
                current = {}
                first_key_of_record = True
                current_key = None
            continue

        if line.startswith("#"):
            continue

        if line.endswith(":"):
            raw_key = line[:-1]

            if first_key_of_record:
                current["name"] = raw_key
                current_key = "atomic_number"
                first_key_of_record = False
            else:
                current_key = raw_key[0]
                if current_key in ("c", "d"):
                    current[current_key] = []
        else:
            values = line.split()
            if current_key == "atomic_number":
                current["atomic_number"] = int(values[0])
            elif current_key in ("c", "d"):
                current[current_key].extend(values)
            else:
                current[current_key] = values[0]

    if current:
        yield current


def atom_stream_ion(lines):
    """begin:/end:-delimited records, with 'Symbol(charge)' notation lines."""
    current = {}
    current_key = None

    for raw_line in lines:
        line = raw_line.strip()

        if not line or line.startswith("#"):
            continue

        if line == "begin:":
            current = {}
            current_key = None
            continue

        if line == "end:":
            if current:
                yield current
            current = {}
            current_key = None
            continue

        if line.endswith(":"):
            label = line[:-1].strip()
            if label == "entry no.":
                current_key = "entry_no"
            elif label.startswith("Z (atomic number"):
                current_key = "atomic_number"
            elif label.startswith("Z0"):
                current_key = "electrons"
            elif label == "Alpha":
                current_key = "a"
            else:
                m = re.match(r'^([a-zA-Z]+)\d+-[a-zA-Z]+\d+$', label)
                base = m.group(1) if m else label
                current_key = base
                if current_key in ("c", "d"):
                    current[current_key] = []
            continue

        # Not a "key:" line -- either a plain value, or the "Symbol(charge)" line
        if current_key == "entry_no":
            if "entry_no" not in current:
                current["entry_no"] = int(line)
            else:
                current["name"] = line
                current_key = None
            continue

        values = line.split()
        if current_key in ("atomic_number", "electrons"):
            current[current_key] = int(values[0])
        elif current_key in ("c", "d"):
            current[current_key].extend(values)
        else:
            current[current_key] = values[0]

    if current:
        yield current







if __name__ == '__main__':
    atoms = []

    neutral_source = "Thorkildsen, G. (2023). Acta Cryst. A79, 318-330; Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Acta Cryst. A79, 59–79."
    with open("neutral_atoms.txt") as f:
        atoms.extend([Atom(d) for d in atom_stream_neutral(f)])
    for atom in atoms:
        atom.source=neutral_source
        atom.upper_limit = 6

    ion_source = "Thorkildsen, G. (2024). Acta Cryst. A80, 129-136; Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Acta Cryst. A79, 229–245."
    with open("ions_1.txt") as f:
        atoms.extend([Atom(d) for d in atom_stream_ion(f)])
    for atom in atoms:
        if not atom.source:
            atom.source=ion_source
            atom.upper_limit = 8



    sorted_atoms = sorted(atoms, key=lambda x: (x.atomic_number, x.electrons))

    with open("atomscat.cif", "w", encoding="utf-8") as file:
        file.write("""#\\#CIF_2.0

data_inv_Mott_Bethe_atom_scattering
loop_
_atom_atype_scat.symbol
_atom_atype_scat.inv_mott_bethe_c
_atom_atype_scat.inv_mott_bethe_as
_atom_atype_scat.inv_mott_bethe_bs
_atom_atype_scat.inv_mott_bethe_lower_limit
_atom_atype_scat.inv_mott_bethe_upper_limit
_atom_atype_scat.source\n""")
        for atom in sorted_atoms:
             file.write(f"{atom.cif_line()}\n")

    # oxygen = sorted_atoms[23]
    # print(oxygen.f_s(0.01))


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
