#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2022-05-16
Purpose: This should be the Tm calculator function as seen in Primer Express 3.0.1 for generating qPCR probes.
"""

import math

class Probe():
    """
    A class to represent TaqManMGB probes.

    Attributes (that you care about):
        sequence (str): The string sequence of the probe
        GC (float): The GC% of the probe
        Tm (float): The Tm (°C) of the probe

    """

    def __init__(self, sequence: str) -> None:
        """
        Construct a TaqManMGB probe.

        Parameters:
            sequence (str): Probe sequence as a string
        """

        self.sequence: str = self.AsciiStringToDNA(sequence)
        self.Tm: float = self.getMGBTmFromStr(self.sequence)
        self.GC: float = (self.sequence.count('G') + self.sequence.count('C'))/len(self.sequence)


    def __repr__(self) -> str:
        return f'probe_seq: {self.sequence}, probe_tm {round(self.Tm, 0)}'


    def CalcTempDop(self, Mg: float, Na: float) -> float:
        """
        Calculates some kind of salt-related condition that alters the Tm of the probe.

        Parameters:
            Mg (float): Magnesium concentration
            Na (float): Sodium concentration

        Returns:
            tempDop (float): Slight alteration to Tm due to salt concentrations.
        """

        dna: float = 0.0
        dmg: float = 0.0

        Tmgdop: dict = {}
        Tnadop: dict = {}

        tempDop: float = 0.0

        Nalog: float = math.log10(Na)
        Mglog: float = math.log10(Mg)

        Tmg: dict = {0: 0.0, 1: 4.0E-5, 2: 2.0E-4, 3: 0.001, 4: 0.005, 5: 0.025}

        Tna: dict = {0: 0.0, 1: 0.02, 2: 0.04, 3: 0.11, 4: 0.21, 5: 0.41, 6: 1.0}

        Tmglog: dict = {0: 0.0, 1: -4.39794, 2: -3.69897, 3: -3.0, 4: -2.30103, 5: -1.60206}

        Tnalog: dict = {0: 0.0, 1: -1.69897, 2: -1.39794, 3: -0.958607, 4: -0.677781, 5: -0.387216, 6: 0.0}

        Ttable = {
            0: {0: -22.9, 1: -22.9, 2: -13.7, 3: -4.9, 4: -0.6, 5: -0.1},
            1: {0: -17.3, 1: -17.3, 2: -11.7, 3: -4.9, 4: -0.3, 5: 0.3},
            2: {0: -11.7, 1: -11.7, 2: -9.7, 3: -4.9, 4: 0.0, 5: 0.7},
            3: {0: -3.2, 1: -3.2, 2: -3.2, 3: -2.4, 4: 0.3, 5: 1.1},
            4: {0: 1.1, 1: 1.1, 2: 1.1, 3: 1.1, 4: 1.3, 5: 1.5},
            5: {0: 4.7, 1: 4.5, 2: 4.3, 3: 4.1, 4: 4.1, 5: 4.1},
            6: {0: 9.6, 1: 9.03, 2: 8.56, 3: 8.1, 4: 7.83, 5: 7.56}
        }

        if (Mg == 0.005 and Na == 0.04):
            tempDop = 0.0
            return tempDop

        img: int = 0
        while img < 5:
            if Mg < Tmg[img + 1]: break
            else: img += 1

        ina: int = 0
        while ina < 6:
            if Na < Tna[ina + 1]: break
            else: ina += 1

        if ina == 6:
            i: int = 0
            while i < 6:
                ad: float = (Ttable[6][i] - Ttable[5][i]) / (Tnalog[6] - Tnalog[5])
                bd: float = Ttable[6][i] - ad * Tnalog[6]
                Tnadop[i] = ad * Nalog + bd
                i += 1

        if img == 5:
            i: int = 0
            while i < 7:
                ad: float = (Ttable[i][5] - Ttable[i][4]) / (Tmglog[5] - Tmglog[4])
                bd: float = Ttable[i][5] - ad * Tmglog[5]
                Tmgdop[i] = ad * Mglog + bd

        if ina == 6 and img == 5:
            ad: float = (Tnadop[5] - Tnadop[4]) / (Tmglog[5] - Tmglog[4])
            bd: float = Tnadop[5] - ad * Tmglog[5]
            tup: float = ad * Mglog + bd
            ad = (Tmgdop[6] - Tmgdop[5]) / (Tnalog[6] - Tnalog[5])
            bd = Tmgdop[6] - ad * Tnalog[6]
            tdown: float = ad * Nalog + bd
            tempDop = (tup + tdown) / 2.0

        if img != 5:
            if img == 0: dmg = Mg / Tmg[1]
            else: dmg = (Mglog - Tmglog[img]) / (Tmglog[img + 1] - Tmglog[img])

        if ina != 6:
            if ina == 0: dna = Na / Tna[1]
            else: dna = (Nalog - Tnalog[ina]) / (Tnalog[ina + 1] - Tnalog[ina])

        if ina != 6 and img != 5:
            tup: float = Ttable[ina][img] + dmg * (Ttable[ina][img + 1] - Ttable[ina][img])
            tdown: float = Ttable[ina + 1][img] + dmg * (Ttable[ina + 1][img + 1] - Ttable[ina + 1][img])
            tempDop = tup + dna * (tdown - tup)

        if img == 5 and ina != 6:
            tempDop = Tmgdop[ina] + dna * (Tmgdop[ina + 1] - Tmgdop[ina])
        if img != 5 and ina == 6:
            tempDop = Tnadop[img] + dmg * (Tnadop[img + 1] - Tnadop[img])

        #print(f'img: {img}, ina: {ina}')
        return tempDop


    def CalcT(self, Mg: float, Na: float, Concen: float, sProba: dict):
        """
        Calculates the Tm (°C) for a given TaqManMGB probe.

        Parameters:
            Mg (float): Magnesium concentration
            Na (float): Sodium concentration
            Concen (float): Nucleic acid concentration
            sProba (dict): Probe string array (it might be set to be a dict right now)

        Returns:
            Tm (float): Raw Tm (°C) for the TaqManMGB probe
        """

        dS = {
            0: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
            1: {0: 0.0, 1: -22.2, 2: -22.4, 3: -21.0, 4: -20.4},
            2: {0: 0.0, 1: -22.7, 2: -19.9, 3: -27.2, 4: -21.0},
            3: {0: 0.0, 1: -22.2, 2: -24.4, 3: -19.9, 4: -22.4},
            4: {0: 0.0, 1: -21.3, 2: -22.2, 3: -22.7, 4: -22.2}
        }

        dH = {
            0: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
            1: {0: 0.0, 1: -7900.0, 2: -8400.0, 3: -7800.0, 4: -7200.0},
            2: {0: 0.0, 1: -8500.0, 2: -8000.0, 3: -10600.0, 4: -7800.0},
            3: {0: 0.0, 1: -8200.0, 2: -9800.0, 3: -8000.0, 4: -8400.0},
            4: {0: 0.0, 1: -7200.0, 2: -8200.0, 3: -8500.0, 4: -7900.0}
        }

        dSPCR = {
            0: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
            1: {0: 0.0, 1: -22.251, 2: -22.652, 3: -17.274, 4: -22.988},
            2: {0: 0.0, 1: -20.343, 2: -20.071, 3: -20.68, 4: -17.274},
            3: {0: 0.0, 1: -23.959, 2: -22.041, 3: -20.071, 4: -22.652},
            4: {0: 0.0, 1: -21.767, 2: -23.959, 3: -20.343, 4: -22.251}
        }

        dHPCR = {
            0: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
            1: {0: 0.0, 1: -7848.149736, 2: -8454.769827, 3: -6389.315578, 4: -8184.491174},
            2: {0: 0.0, 1: -7530.587908, 2: -7858.3913, 3: -8083.611101, 4: -6389.315578},
            3: {0: 0.0, 1: -8715.696056, 2: -8954.301519, 3: -7858.3913, 4: -8454.769827},
            4: {0: 0.0, 1: -7231.831033, 2: -8715.696056, 3: -7530.587908, 4: -7848.149736}
        }

        InitS = {0: 0.0, 1: 4.1, 2: -2.8, 3: -2.8, 4: 4.1}

        InitH = {0: 0.0, 1: 2300.0, 2: 100.0, 3: 100.0, 4: 2300.0}

        InitSPCR = {0: 0.0, 1:-7.0977, 2: -19.002, 3: -19.002, 4: -7.0977}

        InitHPCR = {0: 0.0, 1: -1058.65748, 2: -5431.818926, 3: -5431.818926, 4: -1058.65748}

        MGBPCR = {
            0: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
            1: {0: 0.0, 1: 3.407506222, 2: 1.441885745, 3: 0.664524726, 4: 3.060091908},
            2: {0: 0.0, 1: 0.607476337, 2: 0.90513791, 3: -1.104146165, 4: 2.25315143},
            3: {0: 0.0, 1: 2.542327963, 2: 1.819880026, 3: 1.288047464, 4: 2.004729958},
            4: {0: 0.0, 1: 2.462867141, 2: 2.870499201, 3: -8.4998E-5, 4: 3.3130226}
        }

        t: float = 0.0
        w2: float = 0.0
        wmax: float = 0.0
        lProba: int = len(sProba)
        tempDop: float = 0.0

        if (Na == 0.04 and Mg == 0.005):
            StandardPCR = 1
            NewTable = 1

        else:
            StandardPCR = 0
            NewTable = 0

            if Mg * 100.0 + 1.0E-14 > Na:
                NewTable = 1
                if Mg < 1.0E-14: Mg = 1.0E-14
                if Na < 1.0E-14: Na = 1.0E-14
                tempDop = self.CalcTempDop(Mg, Na)

        if Mg < 1.0E-14: Mg = 1.0E-14
        if Na < 1.0E-14: Na = 1.0E-14

        SS: float = 0.0
        SH: float = 0.0

        i: int = 0
        while i < lProba - 1:
            k: int = sProba[i]
            e2: int = sProba[i + 1]
            if NewTable == 0:
                SS += dS[k][e2]
                SH += dH[k][e2]
            else:
                SS += dSPCR[k][e2]
                SH += dHPCR[k][e2]
            i += 1

        if NewTable == 0:
            SS += InitS[sProba[0]] + InitS[sProba[lProba - 1]]
            SH += InitH[sProba[0]] + InitH[sProba[lProba - 1]]
        else:
            SS += InitSPCR[sProba[0]] + InitSPCR[sProba[lProba - 1]]
            SH += InitHPCR[sProba[0]] + InitHPCR[sProba[lProba - 1]]

        i: int = 0
        j: int = lProba - 6
        dds1: float = 0.0
        while i < 5:
            ddsindex1: int = sProba[j]
            ddsindex2: int = sProba[j + 1]
            ddsvalue: float = MGBPCR[ddsindex1][ddsindex2]
            dds1 += ddsvalue
            i += 1
            j += 1

        e1: int = 0
        if sProba[lProba - 7] == 1: e1 += 100
        if sProba[lProba - 6] == 1: e1 += 10
        if sProba[lProba - 5] == 1: e1 += 1
        #print(sProba[lProba - 7], sProba[lProba - 6], sProba[lProba - 5], e1)
        #print(self.sequence[lProba - 7], self.sequence[lProba - 6], self.sequence[lProba - 5], e1)

        wmax = 0.0
        if e1 == 10: wmax = 1.787172519
        elif e1 == 100: wmax = 1.787172519
        elif e1 == 11: wmax = 1.787172519 + 1.987 * math.log(2.0)
        elif e1 == 101: wmax = 1.787172519 + 1.987 * math.log(2.0)
        elif e1 == 110: wmax = 1.787172519 + 1.987 * math.log(2.0)
        elif e1 == 111: wmax = 1.787172519 + 1.987 * math.log(3.0)

        dds1 += wmax
        if (sProba[lProba - 1] == 2 and sProba[lProba - 2] == 2 and sProba[lProba - 3] == 2):
            dds1 += 3.0

        dds2: float = 0.0
        i: int = 0
        j: int = lProba - 7
        while i < 5:
            dds2 += MGBPCR[sProba[j]][sProba[j + 1]]
            i += 1
            j += 1

        e1 = 0
        if sProba[lProba - 8] == 1: e1 += 100
        if sProba[lProba - 7] == 1: e1 += 10
        if sProba[lProba - 6] == 1: e1 += 1
        #print(sProba[lProba - 8], sProba[lProba - 7], sProba[lProba - 6], e1)
        #print(self.sequence[lProba - 8], self.sequence[lProba - 7], self.sequence[lProba - 6], e1)

        wmax = 0.0
        if e1 == 10: wmax = 1.787172519
        elif e1 == 100: wmax = 1.787172519
        elif e1 == 11: wmax = 1.787172519 + 1.987 * math.log(2.0)
        elif e1 == 101: wmax = 1.787172519 + 1.987 * math.log(2.0)
        elif e1 == 110: wmax = 1.787172519 + 1.987 * math.log(2.0)
        elif e1 == 111: wmax = 1.787172519 + 1.987 * math.log(3.0)

        dds2 += wmax
        if (sProba[lProba - 1] == 2 and sProba[lProba - 2] == 2 and sProba[lProba - 3] == 2):
            dds2 += 2.0

        if dds1 > dds2:
            w2 = dds1
        else:
            w2 = dds2

        if StandardPCR == 1:
            t = SH / (SS + 1.987 * math.log(Concen / 2.0) + w2) - 273.15
        elif NewTable == 0:
            t = SH / (SS + 1.987 * math.log(Concen / 2.0) + w2) - 273.15 + 15.0 * math.log10(Na + 20.0 * Mg)
        else:
            tempDop = self.CalcTempDop(Mg, Na)
            #print(tempDop)
            t = SH / (SS + 1.987 * math.log(Concen / 2.0) + w2) - 273.15 + tempDop
        return t


    def AsciiStringToDNA(self, sequence: str) -> str:
        """
        Restricts string sequenced to allowed nucleotide characters (ACMGRSVTWYHKDBN).

        Parameters:
            sequence (str): Input sequence potentially containing trailing whitespace, ?, or uracil, or disallowed characters

        Returns:
            sequence (str): DNA sequence stripped of whitespace, uracil-replaced,and potentially containing ambiguous nucleotides
        """

        strASCII_DNA = "ACMGRSVTWYHKDBN"
        sequence = sequence.upper().replace('?', 'N')
        sequence = sequence.upper().replace('U', 'T')
        sequence = ''.join([i.upper() for i in sequence if i.upper() in strASCII_DNA])
        return sequence


    def getMGBTmFromStr(self, sequence: str) -> float:
        """
        Function takes a sequence as a string input and should return Tm as a float.

        Parameters:
            sequence (str): Raw probe sequence as a string

        Returns:
            Tm (float): Tm (°C) of the probe
        """
        self.Tm: float = 0.0

        self.resTo1234 = {'A': 1, 'C': 2, 'G': 3, 'T': 4}

        if len(sequence) >= 10:
            self.sProba = {}
            abDNASequence = self.AsciiStringToDNA(sequence)
            try:
                index: int = 0
                while index < len(sequence.upper()):
                    self.sProba[index] = self.resTo1234[abDNASequence[index]]
                    index += 1
                #print(self.sProba)
            except: # BARE EXCEPT -- have fun Michael
                # The unpacked program has an exception for errors with secondary structure, but this may not be applicable here.
                pass
            self.Tm = self.CalcT(0.005, 0.1, 2.0000000000000002E-7, self.sProba)
        return self.Tm


list_of_probes: list = [
    'TCTCAGAATTTGGAGGAGGTA',
    'agcgagacgaactacag',
    'aggcatttcgcccgaat',
    'aggcatttcgcccgaat',
    'acgtgtattcggcttcttaga',
    'actagtacgagcattgc',
    'TCAGAATTTGGAGGAGGTA',
    'cttgctcgaaagagctt',
    'tagcatcatcggctacca'
    ]

for probe in list_of_probes:
    print(Probe(probe).GC, Probe(probe))
