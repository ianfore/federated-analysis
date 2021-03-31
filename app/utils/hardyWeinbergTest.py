
class hardyWeinberg:


    def __init__(self, variants):
        # variants is a dict of the form {"benign": "(13, 323339291, 'A', 'T')": {'AA': 1000, 'Aa': 527, 'aa': 2}, ...,
        #                                   "pathogenic": (13,32333221, 'C', 'G')": {'AA': 37, 'Aa': 399', 'aa': 3}, ...,
        #                                   "vus": ...}
        self.variants = variants

    def chiSquareTest(self, c=0.5, criticalValue=3.84):
        # https://en.wikipedia.org/wiki/Hardy-Weinberg_principle
        # since degrees of freedom = 1,  criticalValue = 3.84
        #for t in ['benign', 'pathogenic', 'vus']:
        n = len(self.variants)
        for v in self.variants:
            # 2. calculate p = (2 x Obs(AA) + Obs(Aa)) / (2 x (Obs(AA) + Obs(Aa) + Obs(aa))
            self.variants[v]['p'] = (2 * self.variants[v]['AA'] + self.variants[v]['Aa'])  \
                               / ( 2 * (self.variants[v]['AA'] + self.variants[v]['Aa'] + self.variants[v]['aa']))

            # 3. calculate q = 1 - p
            self.variants[v]['q'] = 1 - self.variants[v]['p']

            # 4. calculate Exp(AA) = p**2 x n
            expAA = n * self.variants[v]['p'] **2

            # 5. calculate Exp(Aa) = 2 x p * q * n
            expAa = 2 * self.variants[v]['p'] * self.variants[v]['q'] * n

            # 6. calculate Exp(aa) = q**2 x n
            expaa = n * self.variants[v]['q'] **2

            # 7. calculate chi-square = sum[ (O - E)**2 / E ]
            if expAA == 0 or expAa == 0 or expaa == 0:
                self.variants[v]['chisquare'] = 0
            else:
                self.variants[v]['chisquare'] = (1.0/expAA) * (abs(self.variants[v]['AA'] - expAA) - c)**2 + \
                                        (1.0/expAa) * (abs(self.variants[v]['Aa'] - expAa) - c)**2 + \
                                        (1.0/expaa) * (abs(self.variants[v]['aa'] - expaa) - c)**2


        return self.variants