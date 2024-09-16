# Parameter-Extraction-for-2nd-Life-Batteries
Included is code that loads, then processes raw cycling data for 8 pouch cells with LMO/graphite chemistry from a retired EV battery pack. The data is processed for the purpose of extracting parameters characterizing an equivalent circuit model. <br />
Upon generating figures for every cell, it can be observed defining a local minimum or maximum for raw data is ambiguous. In the code, I utilize findpeaks(), however, depending from cell to cell, each instance in which findpeaks() is called, the parameters might need to be tuned due wide variability in the HPPC curves between pouch cells. For some extreme cases, I had to hard code in order to ignore any unwanted points findpeaks() identified as a minimum. 

Here is where the public data set comes from: <br />
X. Cui, M. A. Khan, G. Pozzato, S. Singh, R. Sharma, and S. Onori, “Taking second-life batteries from exhausted to empowered using experiments, data analysis, and health estimation,” 2024
