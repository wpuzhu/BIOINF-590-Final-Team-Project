# A6. The neighbourhood is defined at the right scale

**Statement.** The radius (or k) used for "neighbourhood RNA" matches the scale of the biology: juxtacrine contact (~10–20 µm), short-range paracrine (~30–100 µm), or tissue zones (hundreds of µm).

**If violated.** Too small misses the signal. Too large dilutes it, and at large radii the neighbourhood becomes a proxy for **region identity**, which is non-transferable by construction (it encodes *where* rather than *what*).

**How we test it.** E3 radius sweep with UQ on the peak location (bootstrap the argmax over units). A peak at very large radii, combined with poor E5 region transfer, says the model learned region identity rather than a neighbourhood mechanism.
