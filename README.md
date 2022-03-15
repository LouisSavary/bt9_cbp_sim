# bt9_cbp_sim

```
git clone https://github.com/LouisSavary/bt9_cbp_sim.git
cd bt9_cbp_sim/sim/
make
```
Il faut aussi télécharger les traces de simulation : https://jilp.org/cbp2016/framework.html

et pour exécuter une trace (depuis le dossier sim) :

```
./predictor ../traces/SHORT_MOBILE-54.bt9.trace.gz
```

ou alors 
```
cd ../scripts
./doall.pl # ou doit.pl plus rapide, doeval.pl pour le jeu de traces d'évaluation
```
