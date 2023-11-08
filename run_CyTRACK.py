import numpy as np
import cytrack 

args = cytrack.read_args()
if args.cytrack_help:
	cytrack.help()
elif args.get_template:
	cytrack.get_cytrack_inputs_template()
else:
	cytrack.get_cytrack_main(args.parameterfile)
