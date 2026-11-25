# Asymptotics of the Coefficients of $\prod_{k\ge1}(1+kx^k)^k$ (OEIS A266891)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24521853.svg)](https://doi.org/10.5281/zenodo.24521853)

**Preprint (Zenodo DOI):** [https://doi.org/10.5281/zenodo.24521853](https://doi.org/10.5281/zenodo.24521853)

Let $a(n)=[x^n]\prod_{k\ge1}(1+kx^k)^k$ ([A266891](https://oeis.org/A266891)). Vaclav Kotesovec conjectured in 2018 that $\log a(n)\sim n^{2/3}(2\log(3n)-3)/(4\cdot3^{1/3})$. The manuscript proves

$$\log a(n)=\frac{n^{2/3}\bigl(2\log(3n)-3\bigr)}{4\cdot 3^{1/3}}+O\Bigl(\frac{n^{2/3}}{\log n}\Bigr)$$

and the sharper saddle-point statement $\log a(n)=\min_{s>0}\bigl(sn+\log\prod_k(1+ke^{-sk})^k\bigr)+O(n^{1/3}\sqrt{\log n})$. Exact coefficients up to $n=10^4$ indicate that the saddle-point formula with its Gaussian factor is asymptotically exact; this is stated as a conjecture.

## Contents

- `paper.tex`, `paper.pdf`: the manuscript.
- `verify_a266891.py`: recomputes the coefficients exactly, compares them with the OEIS b-file, and evaluates the saddle-point formula.
- `run_10000.log`: output of `python3 verify_a266891.py 10000` (all 10001 b-file terms match).

## Reproducing

```bash
# standard library only
python3 verify_a266891.py 2000      # downloads b266891.txt once; prints ALL CHECKS PASSED
python3 verify_a266891.py 10000     # full b-file range, slower
```

## Status

Preprint, not yet peer-reviewed. Parts of the derivations and code were developed with AI assistance; every computational claim is backed by the script above.

## License

Manuscript text: CC BY 4.0. Code: MIT License.
