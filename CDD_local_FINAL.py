#!/usr/bin/python
# -*- coding: utf-8 -*-

# usage : python3 cdd_local_FINAL.py

from __future__ import division
from subprocess import Popen
from subprocess import PIPE
import numpy as np
import sys
from sympy import *
import re
from fractions import Fraction
from math import *
from math import fmod


################ Fonctions #########################

####### fonction qui parse le fichier de réaction réversible et celui de la matrice de stoechio
def read_file(rvfile, sfile):
    f = open(rvfile)
    rv_lines = f.readlines()
    f.close()

    f = open(sfile)
    s_lines = f.readlines()
    f.close()

    for rv_el in rv_lines:
        rv_el = rv_el.replace("\n", "")
        rv_el = rv_el.replace("\r", "")
        rv = rv_el.split(' ')

    s = []

    for s_el in s_lines:
        s_el = s_el.replace("\n", "")
        s_el = s_el.replace("\r", "")
        s_el = s_el.replace("\t", " ")
        ss = s_el.split(' ')
        news = []
        for e in ss:
            if e != "":
                news.append(e)
        s.append(news)

    return (rv, s)

############## fonction qui compte le nombre de fonctions réversibles
def nb_rev(rv):
    c = 0
    nb = 0
    l_rev = []
    for i in rv:
        if i == "1":
            nb = nb + 1
            l_rev.append(c)
        c = c + 1
    return (nb, l_rev)

### fonction qui dédouble les réactions reversibles
def dedouble(rev, l_rev, rv, s):
    srev = []
    for l in s:
        f = []
        for e in l:
            f.append(e)
        for i in l_rev:
            if '-' in l[i]:
                t = l[i].replace("-", "")
            else:
                if l[i] != "0":
                    t = '-' + l[i]
                else:
                    t = l[i]
            f.append(t)
        srev.append(f)
    return srev

### fonction qui crée le noyau de la matrice de stochio dédoublée
def noyau(srev):
    from fractions import Fraction
    m = Matrix(srev)
    n = m.nullspace()
    nbcol = len(n)
    noyau = []
    for i in range(len(n[0])):
        noyau.append([])
    for colone in n:
        for i in range(len(colone)):
            e = Fraction(str(colone[i]))
            e = float(e)
            noyau[i].append(e)
    return noyau


### fonction qui récupère fichier nullspace et transforme en fichier cdd
def nullInCDD(noy):

    mat = []
    noyau = []
    position = []
    string = ""

    cptl = 0
    for t in noy:
        nbp = 0
        nbn = 0
        nb0 = 0
        nb1 = 0
        pos = 0
        pos1 = -1
        line = []
        line.append(float("0"))
        l = ""
        for i in t:
            l += " "+ str(i)
            if i == 0:
                nb0 += 1
            elif i == 1:
                nb1 += 1
                nbp += 1
                pos1 = pos
            elif i > 0:
                nbp += 1
            else:
                nbn += 1
            line.append(i)
            pos += 1
        if nbp != 0 and nb0 != (len(line)-1):
            string += "0" + l + "\n"
            mat.append(line)
        if nb1 == 1 and (nb0 == len(line)-2):
            if pos1 not in position:
                noyau.append(cptl)
                position.append(pos1)
        cptl += 1

    r = len(mat)
    m = len(mat[0])

    string = str(r) + " " + str(m) + " integer\n" + string
    return r, m, mat, noyau, string

### fonction qui écrit le fichier CDD
def writeFile(string, logfile):
    fileOut = open('nullspace.ine', 'w')
    fileOut.write("* nullspace.ine\n")
    fileOut.write("H-representation\n")
    fileOut.write("begin\n")
    fileOut.write(string)
    fileOut.write("end")
    fileOut.write("\nlogfile_on ")
    fileOut.write(logfile)
    fileOut.close()



def completion(mat, noyau ,efm):

    efm_comp = []
    n = len(mat)
    for e in efm:
        en = [0]*n
        for i in range(n):
            if i in noyau:
                en[i] = e[noyau.index(i)+1]
            else:
                for j in range(len(noyau)):
                    en[i] += e[1+j]*mat[i][1+j]
        efm_comp.append(en)
    return efm_comp


#### créaction du fichier pour cdd local
def creation_file(mat,e,logfile):
    fileOut = open('local.ine', 'w')
    fileOut.write("* local.ine\n")
    fileOut.write("H-representation\n")
    fileOut.write("begin\n")
    mat_bis = []
    mat_not_used = []
    for m in mat:
        somme = 0.0
        for j in range(1,len(m)):
            somme = somme + m[j]*e[j]
        if somme == 0.0:
            mat_bis.append(m)
        else:
            mat_not_used.append(m)
    re = len(mat_bis)
    met = len(mat_bis[0])
    string = str(re) + " " + str(met) + " " + "integer\n"
    fileOut.write(string)
    for m in mat_bis:
        for mm in m:
            fileOut.write(str(int(mm)))
            fileOut.write(" ")
        fileOut.write("\n")
    fileOut.write("end")
    fileOut.write("\nlogfile_on ")
    fileOut.write(logfile)
    fileOut.close()
    return re, met, mat_not_used


#### fonction qui lance cdd
def cdd(file, time):
    command = "timeout " + str(time) + " ./" + CDD + "/cdd " + file
    proc = Popen(command, shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
    stdout_value = proc.stdout.read() + proc.stderr.read()
    stdout_val = str(stdout_value)
    return stdout_val


#### fonction qui regarde si c'est terminé
def cdd_end(file):
    f = open(file)
    lines = f.readlines()
    f.close()
    return len(lines)!=0

#### fonction qui tri le résultat
def read_res(res,m):
    ress = str(res)
    lines = ress.split("\\n")
    efm = []
    import re

    for l in lines:
        l = l.replace("\n", "")
        l = l.replace("\r", "")
        l = l.replace("\t", " ")
        #regex = "( ( |-){1}[0-9]+[.]?[0-9]*){" + str(m) + "}"
        regex = "^(\*|[a-z]|')"
        p = re.compile(regex)
        if p.match(l)==None:
            l.replace("  "," ")
            t = l.split(" ")
            e = []
            e_ok = True
            for i in t:
                if i!="":
                    try :
                        e.append(float(i))
                    except :
                        e_ok = False
            if e_ok : 
                efm.append(e)
    return efm


### fonction supprime les efm non terminaux
def terminal(efm,mat):
    efm_f = []
    direction = []
    for e in efm:
        final = True
        if(e[0]==1):
            final = False
        else:
            for i in mat:
                somme = 0.0
                for j in range(1,len(i)):
                    somme = somme + i[j]*e[j]
                if somme < i[0]:
                    final = False
        if final:
            efm_f.append(e)
        else:
            if e[0]!=1:
                direction.append(e)
    return efm_f, direction

def shootingray(e, dir, mat):
    efm = []
    for d in dir:
        list_t = []
        for m in mat:
            somme = 0
            sommet = 0
            for i in range(1,len(m)):
                if m[i] != 0:
                    somme += m[i]*e[i]
                    sommet += m[i]*d[i]
            if somme!=0 and sommet!=0:
                t = -somme/sommet
                if t>0:
                    list_t.append(t)   ####### Plus tard : trouver une heuristique meilleure pour garder le tmin
        if len(list_t)!=0:
            t_choose = min(list_t)
            new_e = [0.0]
            for i in range(1,len(e)):
                val = e[i] + t_choose*d[i]
                new_e.append(val)
            if new_e != [0]*len(e):
                efm.append(new_e)
    return efm

def ord_zero(m, r, mm):
    colonne = []
    taille = len(m[0])
    for i in range(taille):
        colonne.append([])
    for l in m:
        for j in range(taille):
            colonne[j].append(l[j])
    col_pos = []
    for c in colonne:
        col_pos.append([abs(x) for x in c])
    tup = tuple(col_pos)
    ind = np.lexsort(tup)
    mat_zero = []
    for i in ind:
        mat_zero.append(m[i])
    string_zero = str(r) + " " + str(mm) + " integer\n"
    for ligne in mat_zero:
        for e in ligne:
            string_zero += str(e) + " "
        string_zero += "\n"
    return mat_zero, string_zero

def tri_degre(deg, efms):
    efm_deg = []
    efm_non_deg = []
    for e in efms:
        cpt0 = 0
        for elmt in e:
            if elmt==0:
                cpt0 += 1
        if cpt0 < deg:
            efm_non_deg.append(e)
        else:
            efm_deg.append(e)
    return efm_deg, efm_non_deg


############ fonction qui creer vecteur bit

def bit(efm):
    b = ""
    if len(efm)<32:
        for i in range(len(efm),32):
            b += "0"
    for e in efm:
        if e == 0:
            b += "0"
        else:
            b += "1"
    val = int(b,2)
    return val


#################### VARIABLE DU PROGRAMME #########################

####### LECTURE FICHIER DE PARAMETRE ########

p = open("PARAMETRE.txt")
p_lines = p.readlines()
p.close()
for line in p_lines:
    line = line.replace("\n", "")
    line = line.replace("\r", "")
    l = line.split(' = ')
    if l[0] == "TIME_ESSAI":
        TIME_ESSAI = int(l[1])
    elif l[0] == "TIME":
        TIME = int(l[1])
    elif l[0] == "TIME_LOCAL":
        TIME_LOCAL = int(l[1])
    elif l[0] == "PATH_CDD":
        CDD = l[1]
    elif l[0] == "HOST":
        host_p = l[1]
    elif l[0] == "USER":
        user_p = l[1]
    elif l[0] == "PASSWORD":
        password_p = l[1]
    elif l[0] == "DATABASE":
        database_p = l[1]
    elif l[0] == "RV_FILE":
        rv_file = l[1]
    elif l[0] == "S_FILE":
        s_file = l[1]
    elif l[0] == "ORDER":
        order_p = l[1]
    elif l[0] == "OUT":
        out_p = l[1]



file = "nullspace.ine"


####################################### PARTIE MYSQL ########################################

try :
    import mysql.connector
except:
    print("Impossible de faire l'import")

conn = mysql.connector.connect(host=host_p,user=user_p,password=password_p, database=database_p)
cursor = conn.cursor()

try :
    cursor.execute("DROP TABLE Efm;")
except:
    pass


####################################### FIN PARTIE MYSQL ########################################



################## Instructions #########################

# lecture des fichiers et extraction des informations utiles
[rv, s] = read_file(rv_file,s_file)
(rev, l_rev) = nb_rev(rv)
srev = dedouble(rev,l_rev,rv,s)
n = noyau(srev)

[r, m, mat, noyau, string] = nullInCDD(n)
[mat_zero, string_zero] = ord_zero(mat,r,m)


# servira a decouper le vecteurs bits en segments de 32 bits
size_long = m-1
size = []
while size_long >= 32:
    size.append(32)
    size_long -= 32
if size_long != 0:
    size.append(size_long)


# creation table sql pour stocker efms

create_table = "CREATE TABLE IF NOT EXISTS Efm ("
insert_table = "INSERT INTO Efm ("
insert_bis = "VALUES ("
select_table = "SELECT * FROM Efm WHERE "
select_table2 = ""

var = "bit"
pk = ""
for i in range(len(size)):
    v = var + str(i+1)
    create_table += v + " int(10) UNSIGNED NOT NULL, "
    insert_table += v + ","
    insert_bis += "%s,"
    select_table += v + "=%s AND "
    pk += v + ","
pk = pk[:-1]
select_table = select_table[:-5]
create_table += "list_pos int(10) DEFAULT NULL, PRIMARY KEY(" + pk + "));"
insert_table += "list_pos) " + insert_bis + "%s);"

cursor.execute(create_table)
cursor.execute("CREATE INDEX IndexEfm USING BTREE ON Efm(bit1);")


# choix de l'ordre des contraintes

if order_p == "order_by_zero":
    ordre_pref = "lexmin"
    ordre2 = 1
elif order_p == "libre":
    ordre = ["mincutoff","maxcutoff","mixcutoff","lexmin","lexmax","minindex"]
    ordre2 = 0
    ordre_pref = ""
    efm_garde = []

    for o in ordre:
        writeFile(string,o)
        stdout_value = cdd(file, TIME_ESSAI)
        efm_prov = read_res(stdout_value,m)
        [efm_t, _] = terminal(efm_prov,mat)
        if len(efm_t) > len(efm_garde):
            efm_garde = efm_t
            ordre_pref = o
    writeFile(string_zero,"lexmin")
    stdout_value = cdd(file, TIME_ESSAI)
    
    efm_prov = read_res(stdout_value,m)
    [efm_t, _] = terminal(efm_prov,mat)
    if len(efm_t) > len(efm_garde):
        efm_garde = efm_t
        ordre_pref = "lexmin"
        ordre2 = 1
else:
    ordre2 = 0
    ordre_pref = order_p


# CDD global :

if ordre2 == 0:
    writeFile(string,ordre_pref)
else:
    writeFile(string_zero,ordre_pref)

stdout_value = cdd(file, TIME)
efm_prov = read_res(stdout_value,m)

# on recupere les graines, on les stocke dans la base de données et dans dictionnaire triées par incidence

[efm_notused2, _] = terminal(efm_prov,mat)
efm_notused = {}
nb_efm_not_used = 0

for i in range(len(efm_notused2)):
    ict = 1 # on exclut la 1er coordonnée de l'efm qui est juste propre à cdd
    tupe = []
    for dec in size:
        b_e = bit(efm_notused2[i][ict:ict+dec])
        ict += dec
        tupe.append(b_e)
    tupe.append(i)
    tupe = tuple(tupe)
    cursor.execute(insert_table,tupe)
    deg_degener = 0
    for elmt in efm_notused2[i]:
        if elmt == 0:
            deg_degener += 1
    if deg_degener not in efm_notused:
        efm_notused[deg_degener] = [efm_notused2[i]]
    else:
        efm_notused[deg_degener].append(efm_notused2[i])
    nb_efm_not_used += 1



nb_cdd = 0
nb_efms = nb_efm_not_used
liste_efms = [] # utile si parametre out = values

fileEnd = open('logfile.txt', 'w')
txt = "Ordre choisi : " + ordre_pref + "\nNombre Efms après CDD global : " + str(nb_efms) + "\nCDD local :\n"
fileEnd.write(txt)


### on fait du local tant qu'on a des graines non visitées
while nb_efm_not_used != 0:
    nb_nvx = 0
    for k in efm_notused.keys():
        if efm_notused[k] != []:
            # enlever un efm de la pile à la fin
            e = efm_notused[k].pop()
            nb_efm_not_used -= 1
            #k_e = k # incidence de la graine en cours
            break
    # il faut mettre à jour la BD pour passer la graine en cours en visité et donc supprimer son indice de efm_notused de la table
    ite = [nb_efm_not_used]
    ite = tuple(ite)
    # l'ajouter aux efm terminaux
    if out_p == "values":
        liste_efms.append(e)

    # creation du fichier cdd local
    if ordre2 == 0:
        [rbis, mbis, mat_not_u] = creation_file(mat,e,ordre_pref)
    else:
        [rbis, mbis, mat_not_u] = creation_file(mat_zero,e,ordre_pref)
    # lancer CDD localement sur cet efm
    stdout_value = cdd("local.ine", TIME_LOCAL)
    nb_cdd += 1
    # récupérer les efms de ce CDD local
    efm_provi = read_res(stdout_value,mbis)
    # trier les terminaux
    [efm_t, dir] = terminal(efm_provi,mat)
    # shooting ray
    efm_t_dir = shootingray(e,dir,mat_not_u)

    efm_t = efm_t + efm_t_dir
    
    # vérifier si dans efm_ter ou dans efm_not used sinon les y ajouter
    for e_l in efm_t:
        ict = 1
        tup = []
        for dec in size:
            b_el = bit(e_l[ict:ict+dec])
            ict += dec
            tup.append(b_el)
        tupe = tuple(tup)
        cursor.execute(select_table,tupe)
        rows = cursor.fetchall()
        if len(rows)==0: # on ajoute l'efm
            tup.append(nb_efm_not_used)
            tupe2 = tuple(tup)
            cursor.execute(insert_table,tupe2)
            deg_degener = 0
            for elmt in e_l:
                if elmt == 0:
                    deg_degener += 1
            if deg_degener not in efm_notused:
                efm_notused[deg_degener] = [e_l]
            else:
                efm_notused[deg_degener].append(e_l)
            nb_efm_not_used += 1
            nb_nvx += 1

    nb_efms += nb_nvx
    txt = str(nb_cdd) + " : " + str(nb_efms) + "\n"
    fileEnd.write(txt)


fileEnd.close()
fileEfms = open('efm.ext', 'w')

if out_p == "values":
    # Completion à partir du noyau
    efm_complet = completion(mat, noyau, liste_efms)
    for e in efm_complet:
        for j in range(0, len(e)):
            fileEfms.write(str(int(e[j])))
            fileEfms.write(" ")
        fileEfms.write("\n")
else:
    cursor.execute("SELECT * FROM Efm;")
    rows = cursor.fetchall()
    for row in rows:
        for elmt in row:
            fileEfms.write(str(elmt))
            fileEfms.write(" ")
        fileEfms.write("\n")

fileEfms.close()
conn.close()