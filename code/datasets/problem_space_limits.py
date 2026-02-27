#! /usr/bin/env python3
"""
File:           problem_space_limits.py
Description:    Create an overview of the number of proteins that can be
                generated for each combination of length and H-ratio.
"""

from math import ceil, factorial


if __name__ == "__main__":
    lengths = [10, 15, 20, 25, 30]

    total_lengths = 0
    total_lengths_generated = 0

    for l in lengths:  # noqa: E741
        l_str = "LENGTH "
        print(f"{l_str}{l}\n" + "=" * ceil((len(l_str) + (l + 0.1) / 10)))
        total = 0
        total_generated = 0

        for i in range(l + 1):
            perms = int(factorial(l) / (factorial(i) * factorial(l - i)))
            total += perms
            generated = min(perms, 1000)
            total_generated += generated
            print(
                f"\tH={str(i).ljust(2)}: "
                + f"({str(round(generated / perms * 100, 4)).ljust(7)}%) {perms}"
            )

        print(f"\n\tTotal:      {total}")
        print(f"\tGenerated:  {total_generated}\n")
        total_lengths += total
        total_lengths_generated += total_generated

    print(f"\nTotal:       {total_lengths}")
    print(f"Generated:   {total_lengths_generated}")
    print(f"Percentage:  {total_lengths_generated / total_lengths * 100}")
