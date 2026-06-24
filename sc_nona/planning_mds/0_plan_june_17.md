Wherever I say "agent", I mean run a separate agent that does research into this specific topic, and create a .md file to report findings to me, but each agent should have full understanding of the context, but also explore on their own terms (I don't want to only see things that fit with my current plan).


I want to incorporate scRNA data into DNA language models, things like NTv3 and Evo2. 

The rationale is: it's hard for models like Evo2 to learn things like promoter-enhancer pairings, being trained on only DNA data. 
- Agent: please look in the literature for benchmarks for promoter-enhancer pairings, and other sorts of tasks that we would like DNA language models to be good at but they're currently bad at. Please produce a CURRENT_LIMS.md for this. Include things like, long-context tasks but also more specific things like the promoter-enhancer. Please make a list of the most important tasks - benchmarks that exist, and tasks that models aren't evaluated on but should be. 

- Agent: please look into the current state of: "are these models good at learning what a "gene" is? and relationships between different genes?" 

The reason that we as humans know things like promoter-enhancer pairs, or about the concept of genes, is... (tell me what you think it is)
- I think one answer is "because we can see things like, RNA transcripts, and causal relationships from perturbation experiments". Therefore, it makes sense to give the model this information. 

I think the purest form of information on the RNA level is scRNA transcripts, full base-pair coverage. 1) should be full base-pair coverage not gene counts (which a) only cover a fraction of the genome b) tell you nothing about splicing (question - does full base-pair coverage tell you something about splicing?)) 2) cannot be 3' or 5' biased. I believe this narrows down our data stream to smart-seq scRNA. (Agent: interrogate this claim). 

And among smart-seq scRNA datasets, (agent: find me datasets that fulfill these properties):
- I want variation over genotype, because that is the point of my DNA language model. 
Other axes that matter are: how many cells in the scRNA samples; how many tissues; how many species. But the important thing is just "variation over genotypes". 

How to use this scRNA data (BigWig style) for the model? 
I'm motivated by the idea of MSAs: passing in a large aligned matrix, or aligned "tracks". And running attention along the length of the sequence, and also vertically in some way. Because a big part of what I want to learn here is "covariation between RNA expression of different genes". 
Agent: investigate other ways of incorporating this information. 