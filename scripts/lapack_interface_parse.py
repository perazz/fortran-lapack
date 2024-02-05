import os
import re

DEBUG = True

def parse_interfaces(file_name):

    # Read whole file
    with open(os.path.join('.',file_name), 'r') as file:
            # Create an empty list to store the lines

            s = []
            d = []
            q = []
            c = []
            z = []
            w = []

            interfaces = []
            all_funs = []

            # Iterate over the lines of the file
            for line in file:

                ls = line.strip().lower()

                # Extract function mame
                m = re.search(r'^public :: stdlib_(.+)$',ls)


                if not m is None:

                    name = m.group(1)

                    if name[0]=='s':
                        s.append(name)
                    elif name[0]=='d':
                        d.append(name)
                    elif name[0]=='q':
                        q.append(name)
                    elif name[0]=='c':
                        c.append(name)
                    elif name[0]=='z':
                        z.append(name)
                    elif name[0]=='w':
                        w.append(name)

                    all_funs.append(name)

                    # Strip initial
                    stripped = name[1:]

                    # Add to interface
                    if not stripped in interfaces:
                        if ('c'+stripped in c and \
                            'z'+stripped in z and \
                            'w'+stripped in w) or \
                           ('s'+stripped in s and \
                            'd'+stripped in d and \
                            'q'+stripped in q):
                           interfaces.append(stripped)

                #if DEBUG: print("raw =" + line.strip())

    for j in range(len(interfaces)):
        occurrence = 0
        for i in range(len(all_funs)):
           if all_funs[i][1:]==interfaces[j]: occurrence += 1
        print(interfaces[j] + " found " + str(occurrence))


    print(interfaces)
    exit(1)





parse_interfaces('lapack_interface_candidates.txt')





