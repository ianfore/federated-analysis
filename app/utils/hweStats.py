import subprocess
import scipy.stats as stats
import hail
import argparse

class hweStats():

    def hailHWTest(AA, Aa, aa):
        return hail.eval(hail.hardy_weinberg_test(AA, Aa, aa)).p_value

    def getFisherExact(a, b, c, d):
        # create 2x2 contingency table

        # run fisher exact test
        oddsRatio, pValue = stats.fisher_exact([[a, b], [c, d]])

        # get confidence interval for odds ratio: CI = e^(ln(OR) +/- [1.96 * sqrt(1/a + 1/b + 1/c + 1/d)])
        '''logOR = log(oddsRatio)
        SE_logOR = sqrt(1/a + 1/b + 1/c + 1/d)
        lowerBound = logOR - 1.96 * SE_logOR
        upperBound = logOR + 1.96 * SE_logOR
    
        CI95_lower = exp(lowerBound)
        CI95_upper = exp(upperBound)'''

        # use R fisher.test to get confidence interval for LOR
        row1 = str('c(' + str(a) + ',' + str(b) + ')')
        row2 = str('c(' + str(c) + ',' + str(d) + ')')
        CIcmd_1 = str('/usr/bin/Rscript -e "fisher.test(rbind(' + row1 + ',' +  row2 + '))"')
        CIcmd_2 = str('| grep -A1 confidence | sed -n "2,2 p" | xargs')
        CIcmd = str(CIcmd_1 + CIcmd_2)
        lb, ub = subprocess.check_output(CIcmd, shell=True).split()

        # cast result = b'0.002439905 19.594803004\n' to floats
        lb = float(lb)
        ub = float(ub)

        # insert OR and p-value into results dict
        roundedOR = round(oddsRatio, 3)
        roundedPvalue =  round(pValue, 3)
        rounded95CI = (round(lb, 2), round(ub, 2))

        #return roundedOR, roundedPvalue, rounded95CI
        return roundedPvalue


    def hardyWeinbergChiSquareTest(AA, Aa, aa, c=0.5):
        # https://en.wikipedia.org/wiki/Hardy-Weinberg_principle
        AA = int(AA)
        Aa = int(Aa)
        aa = int(aa)
        n= AA + Aa + aa

        # 2. calculate p = (2 x Obs(AA) + Obs(Aa)) / (2 x (Obs(AA) + Obs(Aa) + Obs(aa))
        p = (2 * AA + Aa) / ( 2 * (AA + Aa + aa))

        # 3. calculate q = 1 - p
        q = 1 - p

        # 4. calculate Exp(AA) = p**2 x n
        expAA = n * p **2

        # 5. calculate Exp(Aa) = 2 x p * q * n
        expAa = 2 * p * q * n

        # 6. calculate Exp(aa) = q**2 x n
        expaa = n * q **2

        # 7. calculate chi-square = sum[ (O - E)**2 / E ]
        if expAA == 0 or expAa == 0 or expaa == 0:
            chisquare = 0
        else:
            chisquare = (1.0/expAA) * (abs(AA - expAA) - c)**2 + \
                                    (1.0/expAa) * (abs(Aa - expAa) - c)**2 + \
                                    (1.0/expaa) * (abs(aa - expaa) - c)**2

        return chisquare

def main():
    parser = argparse.ArgumentParser(usage="hweStats --AA AA --Aa Aa --aa aa ")
    parser.add_argument("--AA", dest="AA", help="homozygous ref count", default=None)
    parser.add_argument("--Aa", dest="Aa", help="heterozygous count", default=None)
    parser.add_argument("--aa", dest="aa", help="homozygous alt count", default=None)

    options = parser.parse_args()
    AA = int(options.AA)
    Aa = int(options.Aa)
    aa = int(options.aa)

    print('fisher exact p-value: ' + str(hweStats.hailHWTest(AA, Aa, aa)))
    print('chi-square statistic: ' + str(hweStats.hardyWeinbergChiSquareTest(AA, Aa, aa, 0.5)))

if __name__ == "__main__":
    main()