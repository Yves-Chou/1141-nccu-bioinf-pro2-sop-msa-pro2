def calculate_SoP(input_path, score_path, gopen, gextend):
    # Step 1: Alignment
    sequences = []
    cur = ""
    with open(input_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur != "":
                    sequences.append(cur)
                cur = ""
                continue
            else:
                cur += line
        if cur != "":
            sequences.append(cur)

    aln_len = len(sequences[0])
    if any(len(s) != aln_len for s in sequences):
        raise ValueError("All sequences must have the same aligned length.")

    # Step 2: Substitution Matrix
    matrix = {}
    letters = []
    with open(score_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if not letters:
                letters = parts
                for letter in letters:
                    matrix[letter] = {}
                continue
            row_letter = parts[0]
            for col_letter, score_str in zip(letters, parts[1:]):
                try:
                    score = int(score_str)
                except:
                    score = float(score_str)
                matrix[row_letter][col_letter] = score

    # Step 3: SoP Score
    total_score = 0
    num_seqs = len(sequences)
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seqA = sequences[i]
            seqB = sequences[j]
            pair_score = 0
            gapA_run = False
            gapB_run = False
            for pos in range(aln_len):
                a = seqA[pos]
                b = seqB[pos]

                # Both Alignmnets w/ Gap
                if a == '-' and b == '-':
                    if gapA_run or gapB_run:
                        pair_score += gextend
                    else:
                        pair_score += gopen
                    gapA_run = True
                    gapB_run = True
                    continue

                # Alignment A w/ Gap
                if a == '-' and b != '-':
                    if gapA_run:
                        pair_score += gextend
                    else:
                        pair_score += gopen
                        gapA_run = True
                    gapB_run = False
                    continue

                # Alignment A=B w/ Gap
                if b == '-' and a != '-':
                    if gapB_run:
                        pair_score += gextend
                    else:
                        pair_score += gopen
                        gapB_run = True
                    gapA_run = False
                    continue

                if a in matrix and b in matrix[a]:
                    pair_score += matrix[a][b]
                elif b in matrix and a in matrix[b]:
                    pair_score += matrix[b][a]
                gapA_run = False
                gapB_run = False

            total_score += pair_score

    if isinstance(total_score, float) and total_score.is_integer():
        total_score = int(total_score)
    return total_score
