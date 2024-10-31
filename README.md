### EMSim
This Electromagnetic Simulator is a fun little side project of mine. 
My current job rarely allows me to interact with math and physics, so I wanted to find a way to get back to my roots.
The only functionality currently working is the electric field simulation. Look at the notebook at src/examples.ipynb for a demonstration.
NOTE: The code will run slowly if you increase the 'resolution' by decreasing the spacing. You are effectively increasing the number of 
pixels for which the API has to calculate the electric field value in space. This issue with the speed will be addressed later in the [optimize]
branch. For other upcoming work, look in future_work.txt
