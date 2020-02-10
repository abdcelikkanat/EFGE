## Exponential Family Graph Embeddings

##### Description
We emphasize on exponential family distributions to capture rich interaction patterns between nodes in random walk sequences. We introduce the generic exponential family graph embedding model, that generalizes random walk-based network representation learning techniques to exponential family conditional distributions. We study three particular instances of this model, analyzing their properties and showing their relationship to existing unsupervised learning models. Our experimental evaluation on real-world datasets demonstrates that the proposed techniques outperform well-known baseline methods in two downstream machine learning tasks.

##### Compilation

**1.** You can compile the codes by typing the following command:
```
make all
```
##### Learning Representations

**2.** You can learn the representations of nodes by
```
./efge --corpus CORPUS_FILE --output OUTPUT_FILE --method METHOD_NAME{bern,pois,norm}
```
**3.** To see the detailed parameter settings, you can use
```
./efge --help
```

##### References
A. Celikkanat and F. D. Malliaros, [Exponential Family Graph Embeddings](https://arxiv.org/pdf/1911.09007.pdf), The Thirty-Fourth AAAI Conference on Artificial Intelligence (AAAI-20), New York City, New York, 2020.
