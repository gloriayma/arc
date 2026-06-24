input 1: genotype 

1. (L, )
embedding

2. (L, 16) or whatever
you do some convolutions (downsampling)

3. (L/128, d_ntv)


input 2: RNA tracks
(L, D) (D is the number of scRNA tracks you have)

1. do some dumb thing to downsample into 128bp bins, like mean pooling (claude suggestion)
(L/128, D)

2. "per-cell encoder" - no longer scalar per cell per genome position. now a vector
(L/128, D, d_ntv) - doesn't have to be d_ntv; can choose a smaller dim size for memory, as long as you upproject later 

3. axial attention blocks (x4, from claude i think)
    run attention over the length direction
    run attention over the column direction
(L/128, D, d_ntv)

4. collapse along the D dimension by learning a single query vector of dim (, d_ntv) (these queries match against the per-cell per-genome-spot, so different cells contribute different amounts depending on their d_ntv state) and apply it at each L/128 position
(L/128, d_ntv)


A couple of options of how to combine these two: 
- simplest is additive bias initialized at zero weight on input 2


you have the transformer doing things
(L/128, d_ntv)

then feed it through the ntv-native deconv tower 
also add skips from the RNA representation, somehow
(L, d_ntv)

then you can do heads to do MLM on each track (either the DNA track, or any one of the scRNA tracks) - this is stil ntv-native? because ntv post-training has prediction of continuous count data? 
this needs to be permutation-equivariant (input track -> output track correspondence) (make it a per-cell thing, shared across cells)

the weird thing is, i'm doing this on only human cells so it would make much more sense to have some inductive bias on the human genome structure. and yet the hope of this is to not be human-dependent. that would be a different model 

