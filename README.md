# chollrup
Matlab MEX commands for updating and downdating Cholesky factors after rank one changes

## Introduction

Details behind the implementation are given in [M. Seeger: Low Rank Updates for the Cholesky Decomposition](https://infoscience.epfl.ch/record/161468/files/cholupdate.pdf). If you use this code for scientific work, please cite this paper and provide a link to the code.

**NOTE**: These routines are not our idea, but have been adapted from **LINPACK** code:

```bibtex
@techreport{
  title={LINPACK User's Guide},
  authors={Dongarra, J. and Moler, C. and Bunch, J. and Stewart, G.},
  institution={SIAM},
  year={1979}
}
```

We use the routines dchud, dchdd, dchex, all written by G. Stewart. We applied some modifications to them, without in fact changing any numerical properties. The changes are described in our report linked above, and mainly concern calling BLAS routines whenever possible, and avoiding negative entries on the factor diagonal.

## Installation

- Make sure to first install the [essential](https://github.com/mseeger/essential) package.
- Compile MEX files, by running make. This produces DLL files.
- Copy all DLL and *.m files somewhere into your Matlab path.

## How to use it

Study the Matlab help and have a look at the test programs: testprog1.m for CHOLUPRK1, CHOLDNRK1, testprogex.m for CHOLUPEXCH. Also, read my technical report for all the details.

### FST conventions (from essential)

The Cholesky factor argument LFACT in CHOLUPRK1, CHOLDNRK1 can be either lower or upper triangular. At present, the lower triangular variant runs faster (this is what we use in our work). You have to pass this argument using the convention from FST, part of the essential package. Instead of simply L (n-by-n), pass

```python
{L, [1 1 n n], 'L '}
```

If you use the upper triangular variant (not recommended), it's

```python
{R, [1 1 n n], 'U '}
```

The fst_overview.txt coming with the Essential package will tell you what
this means. CHOLUPEXCH works with upper triangular factors only at present.

### Why would I want to use this? Give me an example!

It is the core computational primitive in many methods which, roughly speaking, do sequential Bayesian posterior updates for a linear or generalized linear model. With "core" we mean: this is where the dominant part of the computation takes place.

Some examples you may have heard about:
- Tipping's Sparse Bayesian Learning (Relevance Vector Machine), at
  least in the online variant proposed by Tipping and Faul, which is
  what most people are using.
- The [Informative Vector Machine](https://papers.nips.cc/paper/2240-fast-sparse-gaussian-process-methods-the-informative-vector-machine.pdf) or other sparse Gaussian process techniques, given they are implemented properly.
- [Bayesian Inference and Optimal Design in the Sparse Linear Model](http://www.jmlr.org/papers/volume9/seeger08a/seeger08a.pdf)

For example, for sparse Bayesian inference, we need to maintain

```python
L*L' = X'*X + DIAG(PIVEC)
GAMVEC = L\(B0VEC + BVEC)
```

All vectors are size N, L is N-by-N. An update at site I changes PIVEC by DELPI*V, BVEC by DELB*V, where
```python
V=ZEROS(N,1); V(I)=1;
```
DELPI, DELB are scalars. The following code will do the job (assuming
ABS(DELPI) is not too small; the update should not be done then):

```python
  % Have the following vectors ready. Do not allocate them anew in the
  % update code!
  %   CVEC (N-by-1), SVEC (N-by-1), TEMPVEC (N-by-1)
  % GAMVEC must be 1-by-N
  if delpi>0
    tempvec=zeros(n,1); tempvec(i)=sqrt(delpi);
    y=delb/sqrt(delpi);
    if choluprk1({l,[1 1 n n],'L '},tempvec,cvec,svec,tempvec,gamvec,y)~=0
      error('Numerical error');
    end
  else
    % If you do not have L\TEMPVEC:
    tempvec=zeros(n,1); tempvec[i]=sqrt(-delpi);
    y=-delb/sqrt(-delpi);
    if choluprk1({l,[1 1 n n],'L '},tempvec,cvec,svec,tempvec,0,gamvec,y)~=0
      error('Numerical error');
    end
    % You may have computed L\TV already, say in TEMPVEC. Here, TV is all
    % 0, except 1 at I. Then, use the following:
    tempvec=sqrt(-delpi)*tempvec;
    y=-delb/sqrt(-delpi);
    if choluprk1({l,[1 1 n n],'L '},tempvec,cvec,svec,tempvec,1,gamvec,y)~=0
      error('Numerical error');
    end
  end
```
