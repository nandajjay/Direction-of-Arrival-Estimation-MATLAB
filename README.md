# Direction-of-Arrival-Estimation-MATLAB
This is matlab based project for estimating the direction of arrival of sound.
Techniques used
* GSS-PHAT
* MUSIC high-res DOA
* Kalman filter tracking
* Delay and sum Beamforming
* MVDR beamforming

This project simulates a 4-microphone array(used simulation because of lack of instruments) and detects from which direction a sound is coming, tracks the angle, and enhances the audio coming from that direction — similar to smart speakers (Alexa/Google Home), robot audition systems, and beamforming microphones.

#Objective of Project
Simulate microphone array signal for a moving sound source, Estimate DOA using GCC-PHAT & MUSIC, and stabilaze it using kalman filter. Also perform Delay-and-Sum and MVDR beamforming towards the tracked direction. Output- produces visualization and audio outputs also generates a GIF animation showing real-time DOA tracking.



# Generated Outputs
Audio
1. mic1_ref.wav
2. beamformed_ds.wav
3. beamformed_mvdr.wav

Plots
1. doa_vs_time.png
2. doa_polar_all.png
3. doa_rose_all.png
4. final_tracked_doa.png
5. waveform_mic1_vs_beamformed.png

Gif
1. doa_animation.gif

#Applications
This project demonstrates the same concepts used in smart speakers(Alexa, Google Home), Audio conferensing sysmtes, Voice-activated devices, Robotics & drones, Surveillance microphone arrays, etc.

# License
Please feel free to use for learning, and research purposes,

#Note
I have used llm based ai model for programming purpose.
