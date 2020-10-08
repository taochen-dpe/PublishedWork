clear all;
close all;
var_a = [1]; % Partition coefficient of the vehicle
var_b = [1]; % Factor to alter the metabolic rate

for vehic = var_a
    global vehic
        for kr = var_b
            global kr
                aromatic_amine_model()
        end 
end

