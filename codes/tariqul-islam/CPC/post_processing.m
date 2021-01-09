function Y = post_processing(BPM)
        
        for i = 5: length(BPM)-3
            X = BPM(i-4:i+3);
            m = mean(X);
            s = std(X);
            if abs(BPM(i)-m)< s
                BPM(i) = BPM(i);
            else
                BPM(i) = m;
                
            end
            BPM(i) = 0.7*BPM(i)+0.15*BPM(i-1)+0.15*BPM(i+1);
        end
        Y = BPM;

end