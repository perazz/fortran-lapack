# Revert the ugly labelled continue statements created by fprettify

def revert_fprettify_labels(source_folder,file_name):

    import re
    import os

    # Load whole file; split by lines
    with open(os.path.join(source_folder,file_name), 'r') as file:
        # Create an empty list to store the lines
        file_body = []

        # Iterate over the lines of the file
        for line in file:
            file_body.append(line.rstrip())


    file.close()

    for i in range(len(file_body)):
       if 'continue' in file_body[i]:
           old = file_body[i]
           nspaces = len(file_body[i-1])-len(file_body[i-1].lstrip())
           file_body[i] = re.sub(r'(\s*)([0-9]+)(\s*)(continue)(\s*)'," "*nspaces+r'\2 \4',old)
       elif 'call linalg_error_handling(err0,err)' in file_body[i]:
           old = file_body[i]
           nspaces = len(file_body[i-1])-len(file_body[i-1].lstrip())
           file_body[i] = " "*nspaces+old.strip() 


    # Write out
    fid = open(os.path.join(source_folder,file_name), 'w')
    for i in range(len(file_body)):
        fid.write("{}\n".format(file_body[i]))
    fid.close()


# Launch script
revert_fprettify_labels("../src","stdlib_linalg_lapack_s.f90")
revert_fprettify_labels("../src","stdlib_linalg_lapack_d.f90")
revert_fprettify_labels("../src","stdlib_linalg_lapack_q.f90")
revert_fprettify_labels("../src","stdlib_linalg_lapack_c.f90")
revert_fprettify_labels("../src","stdlib_linalg_lapack_z.f90")
revert_fprettify_labels("../src","stdlib_linalg_lapack_w.f90")
revert_fprettify_labels("../src","stdlib_linalg_solve.f90")
revert_fprettify_labels("../src","stdlib_linalg_inverse.f90")
revert_fprettify_labels("../src","stdlib_linalg_eye.f90")
revert_fprettify_labels("../src","stdlib_linalg_least_squares.f90")
revert_fprettify_labels("../src","stdlib_linalg_determinant.f90")




