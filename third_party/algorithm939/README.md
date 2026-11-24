# Original Algorithm 939

Gil, Segura and Temme, ACM TOMS 40(3), Article 20 (2014).
DOI: https://doi.org/10.1145/2591004

Official distribution: https://calgo.acm.org/939.zip
Pinned SHA256: `dbf77461a9c02c603f8eb2bfcfcef6e3d1cfbb7d628b54b18f79c4a7e62af9e4`.
`939.zip` and `upstream/` retain the original bytes, including documentation and test data.
`source.json` records the archive and individual file hashes.

The CALGO index states that use is subject to the ACM Software Copyright and
License Agreement: https://www.acm.org/publications/policies/software-copyright-notice .
No project license is applied to these third-party files. The archive includes
its original algorithm attribution and documentation; it does not include a
separate permissive LICENSE file. Retain the ACM notice and source attribution.

The original interface is `marcum(mu,x,y,p,q,ierr)` with
`mu=nu, x=a*a/2, y=b*b/2` and supported domain
`1<=mu<=10000, 0<=x,y<=10000`.
Its accuracy constants remain unchanged. The upstream README states an aimed
relative accuracy close to 1e-11 over this domain; this is not a tunable rtol
or a certified per-evaluation error bound.

Our separate C adapter validates the domain before calling the public routine:
the original routine resets its ierr after its domain checks. The adapter also
handles the exact y=0 boundary. No numerical kernels inside upstream were edited.
