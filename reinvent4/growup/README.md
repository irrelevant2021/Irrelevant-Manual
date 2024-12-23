## Molecules grow from fragments by adjusting the score of molecular volume in RL
  
[[stage.scoring.component]]  
[stage.scoring.component.MolVolume]  
[[stage.scoring.component.MolVolume.endpoint]]  
name = "Moleculer Volume (RDKit)"  
weight = 0.5  ###weight
transform.type = "double_sigmoid"  
transform.high = 450.0  
transform.low = 250.0  ###(450+250)/2=350, the volume of the 350 is the best, and it can be worse for both large and small.  
transform.coef_div = 500.0  
transform.coef_si = 20.0  
transform.coef_se = 20.0  ###the larger the coef_si and coef_se(both sides of 350), the steeper the slope
