# Inverse-Dynamics-Volleyball-Spike
This project simulates a 3-DOF Prosthetic Arm for a Volleyball spike using the Recursive Newton-Euler algorithm in MATLAB and the simulation is provided in 2D Working Model software. A PID controller is added to make the motion of the arm more natural.

## Inverse Dynamics
Inverse dynamics is the method to calculate the forces in a robots joint using a set of joint accelerations. They are crucial for the control and motion planning of a Robot.

## Recursive Newton-Euler ALgorithm
The Recursive Newton-Euler equation works by using the **Theta**, **Theta-dot** and **Theta-double-dot** (i.e. which are the set of joint position, velocities, and acceleration) at a particular joint. The algorithm takes these inputs with some initial inputs to return the **Force** and **Torques** at that particular joint.


## Simulation
The simulation is in Working Model 2D software, where three links and three motor joints are connected. The motor joints input the motor torques.

## PID control
The results after feeding the torques were not perfect due to factors like error in motor, gravity, etc. Hence, I decided to add a PID controller to correct the torques.

### Before PID
![model-video-before-pid](https://user-images.githubusercontent.com/71589098/175041764-d8a7822a-953e-4e29-ac14-38f78f1c2067.gif)

### After PID

![model-video-after-pid](https://user-images.githubusercontent.com/71589098/175041940-366fc6db-b617-43be-b67f-73d558a18873.gif)


