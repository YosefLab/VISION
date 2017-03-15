import os

# Remove existing data/network.json file
os.system("rm data/network.json")

# Load in foundation file
os.system("python cgi-bin/network_json_v3.cgi -f data/org.xlsx")


# Load in subnetworks
os.system("python cgi-bin/network_json_v3.cgi -g data/NMPM\ Network.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/ACE.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/plansource.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/Clear\ the\ Decks.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/Diversity\ and\ Inclusion.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/i2V.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/innovatrium.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/synapse.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/Single\ Decision\ Making.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/Aetion.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/Established\ Products.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/MAP.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/gSpeak.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/Biosimilars.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/gCareer.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/IIS\ Webportal.xlsx")
os.system("python cgi-bin/network_json_v3.cgi -g data/Medical\ Team\ Charter\ Operationalization.xlsx")

# Load in resource file
os.system("python cgi-bin/network_json_v3.cgi -b data/group_resources.xlsx")

# Load in insights file
os.system("python cgi-bin/network_json_v3.cgi -i data/insights.xlsx")

# Load in proccess file
os.system("python cgi-bin/network_json_v3.cgi -p data/nmpm_ideation_process.xlsx")





