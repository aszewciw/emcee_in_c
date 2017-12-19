Some quick, unformatted notes on changes to this version.

There are a lot of small changes with regard to data types, (e.g. using size_t
a lot more for array indexing). Quite frankly, I never know what's actually the
preferred dtype to use.

The more significant changes are these:

1. Added a burn-in phase (though I wrote it as lazily as I possibly could have)
2. Removed the "chain" data structure. In other words, I no longer store the full mcmc.
3. Output data at each step in the chain, which will be useful in case the code
   crashes while running on a cluster.
4. Minor change to the header of the ouput file.