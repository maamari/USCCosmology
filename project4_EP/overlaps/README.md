Requirements
--------------
streamlit, protobuf, numpy, healpy, matplotlib, pandas

Using the tool
--------------
See [video](https://drive.google.com/file/d/1NAmC_RPqxRY_AzblxxcFIF4qD7lW6DLG/view?usp=sharing)<br>

Known issues
--------------
If you run into a toml/decoder error, run:<br>
`cd ~/.streamlit`<br>
`rm config.toml`<br>

<br>
The following contains my work from my time with the USC Theoretical Cosmology Group.
<br><br>
Project 1:<br>
Work conducted regarding the validation and extension of the PSZ2 Galaxy Cluster Catalog. Done at the beginning of my sophomore year. Corresponds to the following publication: https://arxiv.org/abs/1907.06364
<br><br>
Project 2:<br>
Work done to enhance computational efficiency in the Lorentz boosting code, CosmoBoost (https://github.com/syasini/CosmoBoost). This involved exploring numerical instabilities to determine why the code was breaking, implementing a new analytic approach to avoid the computational costs of numerically solving the kernel, and implementing a new chunking approach which allowed for large kernels to be passed into the ODE solver without requiring excessive amounts of RAM. Code ran on the USC supercomputer.
<br><br>
Project 3:<br>
In progress: Fisher matrix construction, likelihood analysis


pip install streamlit
pip install --upgrade protobuf

