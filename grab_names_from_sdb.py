from scilifelab.db.statusdb import ProjectSummaryConnection
pcon = ProjectSummaryConnection(url="tools.scilifelab.se", username="mario", password="MNcv78df!")
project = pcon.get_entry('C.Dixelius_13_01')
for sample in project.get("samples",{}).values():
       print("\t".join([sample.get('scilife_name'),sample.get('customer_name')]))
