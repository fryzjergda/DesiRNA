def check_dot_bracket(ss):
    """
    Checks if the given secondary structure string is in the correct dot-bracket notation.
    It ensures that every opening bracket has a matching closing bracket, and vice versa.

    Parameters:
    ss (str): The secondary structure string in dot-bracket notation.

    Returns:
    list: A list of pairs, where each pair is a list of two indices representing a base pair.

    Raises:
    SystemExit: If the secondary structure string is not in the correct dot-bracket notation.
    """

    db_list = [['(', ')'], ['[', ']'], ['<', '>'], ['{', '}'], ['A', 'a'], ['B', 'b'], ['C', 'c'], ['D', 'd'], ['E', 'e']]
    allowed_characters = '()[]<>{}AaBbCcDdEe.&'

    for i in range(0, len(ss)):
        if ss[i] not in allowed_characters:
            print("Not allowed characters in structures. Check input file.")
            sys.exit()

    stack_list = []
    pairs_list = []

    # stack-pop for all versions of brackets form the db_list
    for i in range(0, len(db_list)):
        for c, s in enumerate(ss):
            if s == db_list[i][0]:
                stack_list.append(c)
            elif s == db_list[i][1]:
                if len(stack_list) == 0:
                    print("There is no opening bracket for nt position " + str(c + 1) + '-' + ss[c])
                    sys.exit()
                elif s == db_list[i][1]:
                    pairs_list.append([stack_list.pop(), c])
        if len(stack_list) > 0:
            err = stack_list.pop()
            print("There is no closing bracket for nt position " + str(err) + '-' + ss[err])
            sys.exit()

    return pairs_list
