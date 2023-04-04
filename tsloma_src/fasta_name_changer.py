import sys

# sys.argv[1]
# >NA19240_chr15_75493250-75514000_19000-21999_out        10
# aaaaattagctttatcttgcaactttcttatctctgttaacaatactagcattcttaaatcattcagtcttaaaagtgaagtcaacactgacttacccatttccttcatctcccagaagcagtcagtcaccaaatcctagagatttttacctttgaaatgtttttggttttattcttttgccctcatcttctccccatcctgttcttcctgctcgatttcactgacacaaagaccacgtcatcttctaagttactgtaacaccctcttagtctactttttgaataggcctttttcccttccaatctactctggaaaccccaccagtttaatctttcttaaaaattacttcatcagaacatagacacaccatgctcactttggcctctattctttcgctcaagctgtctcgtcacctcaatgctattcattttctacccattttcaagtgccacctcttccttgaagctttctctgaccaccgcattccttgtgaagcttcttctcttttgcttttgatggcacagctggatctatacaggttgaccatttaatttatatactcttccatccattattttactttgtaggaattaatgatatctcactaaacatggtcctaaggacaaaggattacaggtgtgagccactgcacctggccaaattctcctttgcacagggctgagcactcagggctgagcactcagcaataactttgaagatacttgtgggctaccaactccactttgatcaagacagtaaccggccgggcgcagggctcacacctgtaatcccaacactttgggaggccgaggcaggtggatcacttgaggtcaggagtccgagaccagcctggccaacatggtaaaacccgtctctactaaaaatatgaaaattagctgggcacgatggcaggtgcctgtaatcccagctactcgggagcctgagttaggagaattgtttgaacccaggaggcggaggttgcattaattgcattaaaccaagatcgcaccattgaactccggcctggggcgacagcaagactccatctcaagaaaaaaaaaaacgtaacaactcaagataattccttatactttctgggtgtttctttcctctctgctttcaggaaagcagggcagtactaatggttaattttg

# sys.argv[2] - cnt_F
# sys.argv[3] - cnt_B
# sys.argv[4] - i番目
# sys.argv[5] - 初回ブロックなら1, それ以降なら2
# sys.argv[6] - 1st consensus なら1, 2nd consensus stepなら2
# sys.argv[7] - 1st consensusなら0, 2nd consensusでhap1なら1, 2nd consensusでhap2なら2

if int(sys.argv[6]) == 1:
	if int(sys.argv[5]) == 1:
		with open(sys.argv[1], mode="r") as f:
			for line in f:
				if line[0] == ">":
					line = line.replace("\n", "")
					col = line.split("\t")
					col_1 = col[0].split("_")	
					name = col_1[0]+"_"+col_1[1]+"_"+col_1[2]+"\t"+sys.argv[2]+","+sys.argv[3]+","+sys.argv[4]+"\t"+col[1]
				else:
					line = line.replace("\n","")
					seq = line
		new = name + "\n" + seq	+ "\n"
		with open(sys.argv[1], mode="w") as f:
			f.write(new)

	else:
		with open(sys.argv[1], mode="r") as f:
			for line in f:
				if line[0] == ">":
					line = line.replace("\n", "")
					col = line.split("\t")
					col_1 = col[0].split("_")
					name = col_1[0]+"_"+col_1[1]+"_"+col_1[2]+"\t"+sys.argv[2]+","+sys.argv[3]+","+sys.argv[4]+"\t"
					cov = ""
					for i in range(2,len(col)):
						cov += col[i] + ","
					cov = cov[:-1]
					name += cov
				else:
					line = line.replace("\n","")
					seq = line
		new = name + "\n" + seq	+ "\n"
		with open(sys.argv[1], mode="w") as f:
			f.write(new)


elif int(sys.argv[6]) == 2:
	# hap1 or hap2
	if int(sys.argv[7]) == 1:
		hap = "hap1"
	elif int(sys.argv[7]) == 2:
		hap = "hap2"

	if int(sys.argv[5]) == 1:
		with open(sys.argv[1], mode="r") as f:
			for line in f:
				if line[0] == ">":
					line = line.replace("\n", "")
					col = line.split("\t")
					col_1 = col[0].split("_")	
					name = col_1[0]+"_"+col_1[1]+"_"+col_1[2]+"_"+hap+"\t"+sys.argv[2]+","+sys.argv[3]+","+sys.argv[4]+"\t"+col[1]
				else:
					line = line.replace("\n","")
					seq = line
		new = name + "\n" + seq	+ "\n"
		with open(sys.argv[1], mode="w") as f:
			f.write(new)

	else:
		with open(sys.argv[1], mode="r") as f:
			for line in f:
				if line[0] == ">":
					line = line.replace("\n", "")
					col = line.split("\t")
					col_1 = col[0].split("_")
					name = col_1[0]+"_"+col_1[1]+"_"+col_1[2]+"_"+hap+"\t"+sys.argv[2]+","+sys.argv[3]+","+sys.argv[4]+"\t"
					cov = ""
					for i in range(2,len(col)):
						cov += col[i] + ","
					cov = cov[:-1]
					name += cov
				else:
					line = line.replace("\n","")
					seq = line
		new = name + "\n" + seq	+ "\n"
		with open(sys.argv[1], mode="w") as f:
			f.write(new)


