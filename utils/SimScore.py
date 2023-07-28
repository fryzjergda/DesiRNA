#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# from truth import truth
import sys
import math;

OPEN_BRACKETS={"(":"0","[":"1","<":"2","{":"3","A":"4","B":"5","C":"6","D":"7","E":"8"}
CLOSE_BRACKETS={")":"0","]":"1",">":"2","}":"3","a":"4","b":"5","c":"6","d":"7","e":"8"}

def pairing_positions(s1):
    final_pairs=[];l_opens=[];top=0;  l_closes=[];
    l_dots = [];
    n = 0;
    for i in range(len(s1)):
        if s1[i] in OPEN_BRACKETS:
            l_opens.append([i,OPEN_BRACKETS[s1[i]]])
        if s1[i] in CLOSE_BRACKETS:
            l_closes.append([i,CLOSE_BRACKETS[s1[i]]])
        if s1[i] == '.' or s1[i] == '-':
            final_pairs.append((i, -1));

    lengh=(len(l_opens)-1)
    for i in range(lengh,-1,-1):
        for j in range(lengh+1):
            if (l_opens[i][-1])==(l_closes[j][-1]) and (l_closes[j][0]>l_opens[i][0]):
                #print((stack[i][0],stack1[j][0]))
                final_pairs.append((l_opens[i][0],l_closes[j][0]))
                final_pairs.append((l_closes[j][0],l_opens[i][0]))
                del(l_opens[i][-1])
                del(l_closes[j][-1])

    return(dict(sorted(final_pairs)));

class SimScore:
    def __init__(self,ref_ss,query_ss):
#        print(ref_ss, "target ss")
#        print(query_ss)
        self.ref_ss = ref_ss
        self.query_ss = query_ss
        # truth.AssertThat(len(self.ref_ss)).Named("size of ref_ss").IsEqualTo(len(self.query_ss));

    def find_basepairs(self):
        self.bp_dict_r = pairing_positions(self.ref_ss);
        self.bp_dict_q = pairing_positions(self.query_ss);
        # truth.AssertThat(len(self.bp_dict_r)).Named("size of bp_list_r").IsEqualTo(len(self.bp_dict_q));
    def cofusion_matrix(self):
        tp = 0;
        fp = 0;
        tn = 0;
        fn = 0;
        for i in range(len(self.bp_dict_r)):
            if self.bp_dict_r[i] == self.bp_dict_q[i] and self.bp_dict_r[i] != -1:
                tp += 1;
            if self.bp_dict_r[i] == self.bp_dict_q[i] and self.bp_dict_r[i] == -1:
                tn += 1;
            if self.bp_dict_r[i] != self.bp_dict_q[i]:
                if self.bp_dict_r[i] == -1:
                    fp += 1;
                else:
                    fn += 1;
        self.conf_mat = (tp, fp, fn, tn);

    def mcc(self):
        tp, fp, fn, tn = self.conf_mat;
        if (tp == 0 and fp ==0 and fn == 0 and tn != 0):
            numerator =1
            denominator = 1
        else:
            numerator = ((tp * tn) - (fp * fn));
            denominator = math.sqrt((tp + fp) *  (tp + fn) * (tn + fn) * (tn + fp));
        epsilon = 0.00001;
        return round(numerator / (denominator + epsilon),3);

    def recall(self):
        tp, fp, fn, tn = self.conf_mat;
        return round((tp/(tp + fp + 0.001)), 3);


    def precision(self):
        tp, fp, fn, tn = self.conf_mat;
        return round((tp/(tp + fn + 0.001)), 3);

    def fscore(self):
        makhraj=self.precision() + self.recall()
        if makhraj <0.001:
            makhraj=0.001

        return round(2*(self.precision() * self.recall()/ (makhraj)),4);

    def mcc_reverse(self):
        return -self.mcc()
    def recall_reverse(self):
        return -self.recall()
    def precision_reverse(self):
        return -self.precision()
    def fscore_reverse(self):
        return -self.fscore()


if __name__ == "__main__":
    ssc = SimScore("(((...)))", "((.(-).))");
    ssc.find_basepairs();
    ssc.cofusion_matrix();
    print(ssc.mcc());
    print(ssc.recall());
    print(ssc.precision());
    print(ssc.fscore());


