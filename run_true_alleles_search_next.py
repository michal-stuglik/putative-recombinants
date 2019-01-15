#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# License: GPL v3

import true_alleles_search_next


def main():

    args = ["-f", "true_all_jMHC.fas",
            "-o", "putative_chimera_trios.txt"]
    
    true_alleles_search_next.main(args)


if __name__ == "__main__":
    main()
