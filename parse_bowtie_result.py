import sys
import re


def main():
    for filename in sys.argv:
        results = []
        for line in open(filename):
            if line.strip().endswith(" overall alignment rate"):
                results.append(float(line.split('%')[0]))

        print filename, results
        print sum(results)/len(results)

if __name__ == "__main__":
    main()

