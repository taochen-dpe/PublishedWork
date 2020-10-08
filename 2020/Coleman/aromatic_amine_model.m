function out = aromatic_amine_model(x,xdata)
import com.comsol.model.*
import com.comsol.model.util.*
% Create the model object on the comsol server
model = ModelUtil.create('Model');
model.modelPath('D:\');
model.label('AHT.mph');
model.comments(['Untitled\n\n']);
model.label('Sample_simulation_AHT');

% Show progress bar
ModelUtil.showProgress(true);

%% Global definitions
% Diffusion coefficients
model.param.set('D_v_S', '(K*T)/(6*pi*eta_c*rs_m_S)', 'Diffusion coefficient of S in vehicle');
model.param.set('D_v_P1', '(K*T)/(6*pi*eta_c*rs_m_P1)[kg/s/m]', 'Diffusion coefficient of P1 in vehicle');
model.param.set('D_v_P2', '(K*T)/(6*pi*eta_c*rs_m_P2)[kg/s/m]', 'Diffusion coefficient of P2 in vehicle');
model.param.set('D_v_P3', '(K*T)/(6*pi*eta_c*rs_m_P3)[kg/s/m]', 'Diffusion coefficient of P3 in vehicle');
model.param.set('D_lm_S', '(2E-9)*(exp(-0.46*(rs_A_S^2)))[m^2/s]', 'Diffusion coefficient S, lipid matrix');
model.param.set('D_lm_P1', '(2E-9)*(exp(-0.46*(rs_A_P1^2)))[m^2/s]', 'Diffusion coefficient P1, lipid matrix');
model.param.set('D_lm_P2', '(2E-9)*(exp(-0.46*(rs_A_P2^2)))[m^2/s]', 'Diffusion coefficient P2, lipid matrix');
model.param.set('D_lm_P3', '(2E-9)*(exp(-0.46*(rs_A_P3^2)))[m^2/s]', 'Diffusion coefficient P3, lipid matrix');
model.param.set('D_c_S', '((exp(-alpha*(S_S^lamda)))/(1+(rs_m_S/(k^0.5))+((rs_m_S^2)/(3*k))))*D_v_S', 'Diffusion coefficient S, corneocytes');
model.param.set('D_c_P1', '((exp(-alpha*(S_P1^lamda)))/(1+(rs_m_P1/(k^0.5))+((rs_m_P1^2)/(3*k))))*D_v_P1', 'Diffusion coefficient P1, corneocytes');
model.param.set('D_c_P2', '((exp(-alpha*(S_P2^lamda)))/(1+(rs_m_P2/(k^0.5))+((rs_m_P2^2)/(3*k))))*D_v_P2', 'Diffusion coefficient P2, corneocytes');
model.param.set('D_c_P3', '((exp(-alpha*(S_P3^lamda)))/(1+(rs_m_P3/(k^0.5))+((rs_m_P3^2)/(3*k))))*D_v_P3', 'Diffusion coefficient P3, corneocytes');
model.param.set('D_sc_S', 'AHT_D_sc_S_fine', 'Diffusion coefficient S, stratum corneum');
model.param.set('D_sc_P1', 'AHT_D_sc_P1_fine', 'Diffusion coefficient P1, stratum corneum');
model.param.set('D_sc_P2', 'AHT_D_sc_P2_fine', 'Diffusion coefficient P2, stratum corneum');
model.param.set('D_sc_P3', 'AHT_D_sc_P3_fine', 'Diffusion coefficient P3, stratum corneum');
model.param.set('D_ve_S', '((10^(-8.15-0.655*(log10(MW_S_nounits))))/(0.68+(0.32/fu_S)+(0.025*fnon_S*K_lm_S)))[(m^2)/s]', 'Diffusion coefficient S, viable epidermis');
model.param.set('D_ve_P1', '((10^(-8.15-0.655*(log10(MW_P1_nounits))))/(0.68+(0.32/fu_P1)+(0.025*fnon_P1*K_lm_P1)))[(m^2)/s]', 'Diffusion coefficient P1, viable epidermis');
model.param.set('D_ve_P2', '((10^(-8.15-0.655*(log10(MW_P2_nounits))))/(0.68+(0.32/fu_P2)+(0.025*fnon_P2*K_lm_P2)))[(m^2)/s]', 'Diffusion coefficient P2, viable epidermis');
model.param.set('D_ve_P3', '((10^(-8.15-0.655*(log10(MW_P3_nounits))))/(0.68+(0.32/fu_P3)+(0.025*fnon_P3*K_lm_P3)))[(m^2)/s]', 'Diffusion coefficient P3, viable epidermis');
model.param.set('D_d_S', '((10^(-8.15-0.655*(log10(MW_S_nounits))))/(0.68+(0.32/fu_S)+(0.025*fnon_S*K_lm_S)))[(m^2)/s]', 'Diffusion coefficient S, dermis');
model.param.set('D_d_P1', '((10^(-8.15-0.655*(log10(MW_P1_nounits))))/(0.68+(0.32/fu_P1)+(0.025*fnon_P1*K_lm_P1)))[(m^2)/s]', 'Diffusion coefficient P1, dermis');
model.param.set('D_d_P2', '((10^(-8.15-0.655*(log10(MW_P2_nounits))))/(0.68+(0.32/fu_P2)+(0.025*fnon_P2*K_lm_P2)))[(m^2)/s]', 'Diffusion coefficient P2, dermis');
model.param.set('D_d_P3', '((10^(-8.15-0.655*(log10(MW_P3_nounits))))/(0.68+(0.32/fu_P3)+(0.025*fnon_P3*K_lm_P3)))[(m^2)/s]', 'Diffusion coefficient P3, dermis');
model.param.set('D_r_S', '(K*T)/(6*pi*eta*rs_m_S)[kg/s/m]', 'Diffusion coefficient of S in receptor');
model.param.set('D_r_P1', '(K*T)/(6*pi*eta*rs_m_P1)[kg/s/m]', 'Diffusion coefficient of P1 in receptor');
model.param.set('D_r_P2', '(K*T)/(6*pi*eta*rs_m_P2)[kg/s/m]', 'Diffusion coefficient of P2 in receptor');
model.param.set('D_r_P3', '(K*T)/(6*pi*eta*rs_m_P3)[kg/s/m]', 'Diffusion coefficient of P3 in receptor');

% Partition coefficients
%model.param.set('K_v_S', '1', 'Partition coefficient of S, vehicle/water, unity');
model.param.set('K_v_P1', '1', 'Partition coefficient of P1, vehicle/water, unity');
model.param.set('K_v_P2', '1', 'Partition coefficient of P2, vehicle/water, unity');
model.param.set('K_v_P3', '1', 'Partition coefficient of P3, vehicle/water, unity');
model.param.set('K_lm_S', '(rho_l/rho_w)*(Kow_S^0.69)', 'Partition coefficient, lipid matrix S');
model.param.set('K_lm_P1', '(rho_l/rho_w)*(Kow_P1^0.69)', 'Partition coefficient, lipid matrix P1');
model.param.set('K_lm_P2', '(rho_l/rho_w)*(Kow_P2^0.69)', 'Partition coefficient, lipid matrix P2');
model.param.set('K_lm_P3', '(rho_l/rho_w)*(Kow_P3^0.69)', 'Partition coefficient, lipid matrix P3');
model.param.set('K_c_S', '((1-phib)*K_k_S)+thetab', 'Partition coefficient S, corneocytes');
model.param.set('K_c_P1', '((1-phib)*K_k_P1)+thetab', 'Partition coefficient P1, corneocytes');
model.param.set('K_c_P2', '((1-phib)*K_k_P2)+thetab', 'Partition coefficient P2, corneocytes');
model.param.set('K_c_P3', '((1-phib)*K_k_P3)+thetab', 'Partition coefficient P3, corneocytes');
model.param.set('K_k_S', '(rho_k/rho_w)*4.2*(Kow_S^0.31)', 'Partition coefficient S, keratin');
model.param.set('K_k_P1', '(rho_k/rho_w)*4.2*(Kow_P1^0.31)', 'Partition coefficient P1, keratin');
model.param.set('K_k_P2', '(rho_k/rho_w)*4.2*(Kow_P2^0.31)', 'Partition coefficient P2, keratin');
model.param.set('K_k_P3', '(rho_k/rho_w)*4.2*(Kow_P3^0.31)', 'Partition coefficient P3, keratin');
model.param.set('K_sc_S', '((0.1867*1.37*4.23*(Kow_S^0.31))+(0.0316*0.9*(Kow_S^0.69))+0.7817)', 'Partition coefficient S, stratum corneum');
model.param.set('K_sc_P1', '((0.1867*1.37*4.23*(Kow_P1^0.31))+(0.0316*0.9*(Kow_P1^0.69))+0.7817)', 'Partition coefficient P1, stratum corneum');
model.param.set('K_sc_P2', '((0.1867*1.37*4.23*(Kow_P2^0.31))+(0.0316*0.9*(Kow_P2^0.69))+0.7817)', 'Partition coefficient P2, stratum corneum');
model.param.set('K_sc_P3', '((0.1867*1.37*4.23*(Kow_P3^0.31))+(0.0316*0.9*(Kow_P3^0.69))+0.7817)', 'Partition coefficient P3, stratum corneum');
model.param.set('K_ve_S', '0.7*(0.68+(0.32/fu_S)+(0.025*fnon_S*Kow_S^0.7))', 'Partition coefficient S, viable epidermis');
model.param.set('K_ve_P1', '0.7*(0.68+(0.32/fu_P1)+(0.025*fnon_P1*Kow_P1^0.7))', 'Partition coefficient P1, viable epidermis');
model.param.set('K_ve_P2', '0.7*(0.68+(0.32/fu_P2)+(0.025*fnon_P2*Kow_P2^0.7))', 'Partition coefficient P2, viable epidermis');
model.param.set('K_ve_P3', '0.7*(0.68+(0.32/fu_P3)+(0.025*fnon_P3*Kow_P3^0.7))', 'Partition coefficient P3, viable epidermis');
model.param.set('K_d_S', '0.7*(0.68+(0.32/fu_S)+(0.025*fnon_S*Kow_S^0.7))', 'Partition coefficient S, dermis');
model.param.set('K_d_P1', '0.7*(0.68+(0.32/fu_P1)+(0.025*fnon_P1*Kow_P1^0.7))', 'Partition coefficient P1, dermis');
model.param.set('K_d_P2', '0.7*(0.68+(0.32/fu_P2)+(0.025*fnon_P2*Kow_P2^0.7))', 'Partition coefficient P2, dermis');
model.param.set('K_d_P3', '0.7*(0.68+(0.32/fu_P3)+(0.025*fnon_P3*Kow_P3^0.7))', 'Partition coefficient P3, dermis');
model.param.set('K_r_S', '1', 'Partition coefficient of S, receptot/water');
model.param.set('K_r_P1', '1', 'Partition coefficient of P1, receptor/water');
model.param.set('K_r_P2', '1', 'Partition coefficient of P2, receptor/water');
model.param.set('K_r_P3', '1', 'Partition coefficient of P3, receptor/water');

% Compound data
model.param.set('conc', '((1.5[mg]/(10^-4[m^2]))/MW_S)/(10^-4[m])', 'initial concentration');
%model.param.set('k_rate', 'AHT_k_rate', 'rate of reaction, 3.8055E-05');
model.param.set('MW_S', 'AHT_MW_S', 'molar mass of substrate');
model.param.set('MW_P1', 'AHT_MW_P1', 'molar mass of product 1');
model.param.set('MW_P2', 'AHT_MW_P2', 'molar mass of product 2');
model.param.set('MW_P3', 'AHT_MW_P3', 'molar mass of product 3');
model.param.set('Kow_S', 'AHT_Kow_S', 'log partition coefficient of S in octanol/water');
model.param.set('Kow_P1', 'AHT_Kow_P1', 'log partition coefficient of P1 in octanol/water');
model.param.set('Kow_P2', 'AHT_Kow_P2', 'log partition coefficient of P2 in octanol/water');
model.param.set('Kow_P3', 'AHT_Kow_P3', 'log partition coefficient of P3 in octanol/water');
model.param.set('fu_S', 'AHT_fu_S', 'fraction of unbound solute S');
model.param.set('fu_P1', 'AHT_fu_P1', 'fraction of unbound solute P1');
model.param.set('fu_P2', 'AHT_fu_P2', 'fraction of unbound solute P2');
model.param.set('fu_P3', 'AHT_fu_P3', 'fraction of unbound solute P3');
model.param.set('fnon_S', 'AHT_fnon_S', 'fraction of un-ionised solute S');
model.param.set('fnon_P1', 'AHT_fnon_P1', 'fraction of un-ionised solute P1');
model.param.set('fnon_P2', 'AHT_fnon_P2', 'fraction of un-ionised solute P2');
model.param.set('fnon_P3', 'AHT_fnon_P3', 'fraction of un-ionised solute P3');
model.param.set('MW_S_nounits', '(MW_S/nounits)*1000');
model.param.set('MW_P1_nounits', '(MW_P1/nounits)*1000');
model.param.set('MW_P2_nounits', '(MW_P2/nounits)*1000');
model.param.set('MW_P3_nounits', '(MW_P3/nounits)*1000');

% Other parameters
model.param.set('T', '305.15[K]', 'Temperature');
model.param.set('R', '8.3144598[J/mol/K]', 'Gas Constant');
model.param.set('K', '1.3806E-23[(kg*m^2)/(K*s^2)]', 'Boltzmann''s constant');
model.param.set('rho_k', '1370[kg/m^3]', 'Bulk density of keratin');
model.param.set('rho_l', '900[kg/m^3]', 'Bulk density of lipid');
model.param.set('rho_w', '1000[kg/m^3]', 'Bulk density of water');
model.param.set('H', '0.319[Pa*m^3/mol]', 'Henry''s Law Constant');
model.param.set('N', '6.022141E23[1/mol]', 'Avogadro''s Number');
model.param.set('eta', '0.000765[kg/m/s]', 'Water dynamic viscosity at 36oC, 0.00071');
model.param.set('eta_c', '5[kg/m/s]', 'Cream dynamic viscosity');
model.param.set('rs_m_S', '(((3/4/pi)*0.9087*MW_S_nounits)^(1/3))*(1E-10)[m]', 'Radius of solute S in metres');
model.param.set('rs_m_P1', '(((3/4/pi)*0.9087*MW_P1_nounits)^(1/3))*(1E-10)[m]', 'Radius of solute P1 in metres');
model.param.set('rs_m_P2', '(((3/4/pi)*0.9087*MW_P2_nounits)^(1/3))*(1E-10)[m]', 'Radius of solute P2 in metres');
model.param.set('rs_m_P3', '(((3/4/pi)*0.9087*MW_P3_nounits)^(1/3))*(1E-10)[m]', 'Radius of solute P3 in metres');
model.param.set('rs_A_S', '(((3/4/pi)*0.9087*MW_S_nounits)^(1/3))', 'Radius of solute S in Angstroms');
model.param.set('rs_A_P1', '(((3/4/pi)*0.9087*MW_P1_nounits)^(1/3))', 'Radius of solute P1 in Angstroms');
model.param.set('rs_A_P2', '(((3/4/pi)*0.9087*MW_P2_nounits)^(1/3))', 'Radius of solute P2 in Angstroms');
model.param.set('rs_A_P3', '(((3/4/pi)*0.9087*MW_P3_nounits)^(1/3))', 'Radius of solute P3 in Angstroms');
model.param.set('phib', '0.650028', 'Volume fraction of water in corneocytes, saturation');
model.param.set('thetab', '0.6520', 'Volume fraction of water in corneocytes, actual');
model.param.set('rf', '35E-10[m]', 'Radius of keratin microfibril');
model.param.set('k', 'beta*(rf^2)*((1-thetab)^gamma)');
model.param.set('S_S', '(1-thetab)*(((rs_m_S+rf)/rf)^2)');
model.param.set('S_P1', '(1-thetab)*(((rs_m_P3+rf)/rf)^2)');
model.param.set('S_P2', '(1-thetab)*(((rs_m_P2+rf)/rf)^2)');
model.param.set('S_P3', '(1-thetab)*(((rs_m_P3+rf)/rf)^2)');

% Correction factors/Constants
model.param.set('lamda', '1.09');
model.param.set('gamma', '-1.17');
model.param.set('alpha', '9.47');
model.param.set('beta', '9.32*10^-8');
model.param.set('Structure', 'parameters');
model.param.set('VF_w', 'phib', 'Volume fraction of water in corneocytes');
model.param.set('VF_k', '1-phib', 'Volume fraction of keratin in corneocytes');
model.param.set('VF_lip_SC', '(70.509375/806.509375)', 'Volume fraction of lipid matrix in stratum corneum');
model.param.set('VF_cor_SC', '(736/806.509375)', 'Volume fraction of corneocytes in stratum corneum');
model.param.set('lambda_S', 'rs_m_S/rf', 'Ratio of hydrodynamic radius of solute to radius of keratin microfibril S');
model.param.set('lambda_P1', 'rs_m_P1/rf', 'Ratio of hydrodynamic radius of solute to radius of keratin microfibril P1');
model.param.set('lambda_P2', 'rs_m_P2/rf', 'Ratio of hydrodynamic radius of solute to radius of keratin microfibril P2');
model.param.set('lambda_P3', 'rs_m_P3/rf', 'Ratio of hydrodynamic radius of solute to radius of keratin microfibril P3');
model.param.set('Phi_p_S', 'VF_k*(1+lambda_S)^2', 'Volume fraction of corneocyte inaccesible to solute S');
model.param.set('Phi_p_P1', 'VF_k*(1+lambda_P1)^2', 'Volume fraction of corneocyte inaccesible to solute P1');
model.param.set('Phi_p_P2', 'VF_k*(1+lambda_P2)^2', 'Volume fraction of corneocyte inaccesible to solute P2');
model.param.set('Phi_p_P3', 'VF_k*(1+lambda_P3)^2', 'Volume fraction of corneocyte inaccesible to solute P3');
model.param.set('nounits', '1[kg/mol]');

% Partition coefficient ratio input for boundaries
model.param.set('K_vsc_S', 'K_v_S/K_sc_S');
model.param.set('K_scv_S', 'K_sc_S/K_v_S');
model.param.set('K_vsc_P1', 'K_v_P1/K_sc_P1');
model.param.set('K_scv_P1', 'K_sc_P1/K_v_P1');
model.param.set('K_vsc_P2', 'K_v_P2/K_sc_P2');
model.param.set('K_scv_P2', 'K_sc_P2/K_v_P2');
model.param.set('K_vsc_P3', 'K_v_P3/K_sc_P3');
model.param.set('K_scv_P3', 'K_sc_P3/K_v_P3');
model.param.set('K_scve_S', 'K_sc_S/K_ve_S');
model.param.set('K_vesc_S', 'K_ve_S/K_sc_S');
model.param.set('K_scve_P1', 'K_sc_P1/K_ve_P1');
model.param.set('K_vesc_P1', 'K_ve_P1/K_sc_P1');
model.param.set('K_scve_P2', 'K_sc_P2/K_ve_P2');
model.param.set('K_vesc_P2', 'K_ve_P2/K_sc_P2');
model.param.set('K_scve_P3', 'K_sc_P3/K_ve_P3');
model.param.set('K_vesc_P3', 'K_ve_P3/K_sc_P3');
model.param.set('K_ved_S', 'K_ve_S/K_d_S');
model.param.set('K_dve_S', 'K_d_S/K_ve_S');
model.param.set('K_ved_P1', 'K_ve_P1/K_d_P1');
model.param.set('K_dve_P1', 'K_d_P1/K_ve_P1');
model.param.set('K_ved_P2', 'K_ve_P2/K_d_P2');
model.param.set('K_dve_P2', 'K_d_P2/K_ve_P2');
model.param.set('K_ved_P3', 'K_ve_P3/K_d_P3');
model.param.set('K_dve_P3', 'K_d_P3/K_ve_P3');
model.param.set('K_dr_S', 'K_d_S/K_r_S');
model.param.set('K_rd_S', 'K_r_S/K_d_S');
model.param.set('K_dr_P1', 'K_d_P1/K_r_P1');
model.param.set('K_rd_P1', 'K_r_P1/K_d_P1');
model.param.set('K_dr_P2', 'K_d_P2/K_r_P2');
model.param.set('K_rd_P2', 'K_r_P2/K_d_P2');
model.param.set('K_dr_P3', 'K_d_P3/K_r_P3');
model.param.set('K_rd_P3', 'K_r_P3/K_d_P3');
model.param.set('t_step', '60[s]');

%% Aromatic amine data
% AHT data
model.param.set('AHT_k_rate', '3.8055E-05[1/s]', 'rate of reaction');
model.param.set('AHT_MW_S', '123.155 [g/mol]', 'molar mass of substrate');
model.param.set('AHT_MW_P1', '165.192 [g/mol]', 'molar mass of product 1');
model.param.set('AHT_MW_P2', '203.216 [g/mol]', 'molar mass of product 2');
model.param.set('AHT_MW_P3', '313.306 [g/mol]', 'molar mass of product 3');
model.param.set('AHT_Kow_S', '10^0.80', 'log partition coefficient of S in octanol/water');
model.param.set('AHT_Kow_P1', '10^1.19', 'log partition coefficient of P1 in octanol/water');
model.param.set('AHT_Kow_P2', '10^0.46', 'log partition coefficient of P2 in octanol/water');
model.param.set('AHT_Kow_P3', '10^-0.09', 'log partition coefficient of P3 in octanol/water');
model.param.set('AHT_fu_S', '0.6588', 'fraction of unbound solute S');
model.param.set('AHT_fu_P1', '0.8334', 'fraction of unbound solute P1');
model.param.set('AHT_fu_P2', '0.0109', 'fraction of unbound solute P2');
model.param.set('AHT_fu_P3', '0.441', 'fraction of unbound solute P3');
model.param.set('AHT_fnon_S', '0.9963', 'fraction of un-ionised solute S');
model.param.set('AHT_fnon_P1', '0.9901', 'fraction of un-ionised solute P1');
model.param.set('AHT_fnon_P2', '0.0000', 'fraction of un-ionised solute P2');
model.param.set('AHT_fnon_P3', '0.0010', 'fraction of un-ionised solute P3');
% AMC data
model.param.set('AMC_k_rate', '2.83333E-05[1/s]', 'rate of reaction');
model.param.set('AMC_MW_S', '123.155 [g/mol]', 'molar mass of substrate');
model.param.set('AMC_MW_P1', '165.192 [g/mol]', 'molar mass of product 1');
model.param.set('AMC_MW_P2', '203.216 [g/mol]', 'molar mass of product 2');
model.param.set('AMC_MW_P3', '313.306 [g/mol]', 'molar mass of product 3');
model.param.set('AMC_Kow_S', '10^0.17', 'log partition coefficient of S in octanol/water');
model.param.set('AMC_Kow_P1', '10^0.8', 'log partition coefficient of P1 in octanol/water');
model.param.set('AMC_Kow_P2', '10^0.18', 'log partition coefficient of P2 in octanol/water');
model.param.set('AMC_Kow_P3', '10^-0.18', 'log partition coefficient of P3 in octanol/water');
model.param.set('AMC_fu_S', '0.8211', 'fraction of unbound solute S');
model.param.set('AMC_fu_P1', '0.8776', 'fraction of unbound solute P1');
model.param.set('AMC_fu_P2', '0.0131', 'fraction of unbound solute P2');
model.param.set('AMC_fu_P3', '0.4768', 'fraction of unbound solute P3');
model.param.set('AMC_fnon_S', '0.989', 'fraction of un-ionised solute S');
model.param.set('AMC_fnon_P1', '0.9984', 'fraction of un-ionised solute P1');
model.param.set('AMC_fnon_P2', '0.0', 'fraction of un-ionised solute P2');
model.param.set('AMC_fnon_P3', '0.0005', 'fraction of un-ionised solute P3');
% AEP data
model.param.set('AEP_k_rate', '2.26111E-04[1/s]', 'rate of reaction');
model.param.set('AEP_MW_S', '137.182 [g/mol]', 'molar mass of substrate');
model.param.set('AEP_MW_P1', '179.219 [g/mol]', 'molar mass of product 1');
model.param.set('AEP_MW_P2', '217.239 [g/mol]', 'molar mass of product 2');
model.param.set('AEP_MW_P3', '327.333 [g/mol]', 'molar mass of product 3');
model.param.set('AEP_Kow_S', '10^1.43', 'log partition coefficient of S in octanol/water');
model.param.set('AEP_Kow_P1', '10^1.71', 'log partition coefficient of P1 in octanol/water');
model.param.set('AEP_Kow_P2', '10^1.08', 'log partition coefficient of P2 in octanol/water');
model.param.set('AEP_Kow_P3', '10^0.7', 'log partition coefficient of P3 in octanol/water');
model.param.set('AEP_fu_S', '0.5464', 'fraction of unbound solute S');
model.param.set('AEP_fu_P1', '0.4908', 'fraction of unbound solute P1');
model.param.set('AEP_fu_P2', '0.0104', 'fraction of unbound solute P2');
model.param.set('AEP_fu_P3', '0.3459', 'fraction of unbound solute P3');
model.param.set('AEP_fnon_S', '0.9967', 'fraction of un-ionised solute S');
model.param.set('AEP_fnon_P1', '0.9921', 'fraction of un-ionised solute P1');
model.param.set('AEP_fnon_P2', '0.0', 'fraction of un-ionised solute P2');
model.param.set('AEP_fnon_P3', '0.0002', 'fraction of un-ionised solute P3');
% TDA data
model.param.set('TDA_k_rate', '1.33333E-05[1/s]', 'rate of reaction');
model.param.set('TDA_MW_S', '122.171 [g/mol]', 'molar mass of substrate');
model.param.set('TDA_MW_P1', '164.208 [g/mol]', 'molar mass of product 1');
model.param.set('TDA_MW_P2', '164.208 [g/mol]', 'molar mass of product 2');
model.param.set('TDA_MW_P3', '206.245 [g/mol]', 'molar mass of product 3');
model.param.set('TDA_Kow_S', '10^-0.22', 'log partition coefficient of S in octanol/water');
model.param.set('TDA_Kow_P1', '10^0.54', 'log partition coefficient of P1 in octanol/water');
model.param.set('TDA_Kow_P2', '10^0.54', 'log partition coefficient of P2 in octanol/water');
model.param.set('TDA_Kow_P3', '10^0.01', 'log partition coefficient of P3 in octanol/water');
model.param.set('TDA_fu_S', '0.6626', 'fraction of unbound solute S');
model.param.set('TDA_fu_P1', '0.7684', 'fraction of unbound solute P1');
model.param.set('TDA_fu_P2', '0.7684', 'fraction of unbound solute P2');
model.param.set('TDA_fu_P3', '0.7753', 'fraction of unbound solute P3');
model.param.set('TDA_fnon_S', '0.9264', 'fraction of un-ionised solute S');
model.param.set('TDA_fnon_P1', '0.998', 'fraction of un-ionised solute P1');
model.param.set('TDA_fnon_P2', '0.998', 'fraction of un-ionised solute P2');
model.param.set('TDA_fnon_P3', '0.1', 'fraction of un-ionised solute P3');
% HDAP data
model.param.set('HDAP_k_rate', '5.33056E-04[1/s]', 'rate of reaction');
model.param.set('HDAP_MW_S', '142.162 [g/mol]', 'molar mass of substrate');
model.param.set('HDAP_MW_P1', '184.199 [g/mol]', 'molar mass of product');
model.param.set('HDAP_Kow_S', '10^-1.52', 'log partition coefficient of S in octanol/water');
model.param.set('HDAP_Kow_P1', '10^-1.52', 'log partition coefficient of P in octanol/water');
model.param.set('HDAP_fu_S', '0.8233', 'fraction of unbound solute S');
model.param.set('HDAP_fu_P1', '0.8233', 'fraction of unbound solute P');
model.param.set('HDAP_fnon_S', '0.9999', 'fraction of un-ionised solute S');
model.param.set('HDAP_fnon_P1', '0.9999', 'fraction of un-ionised solute P');
% PPD data
model.param.set('PPD_k_rate', '1.55556E-05[1/s]', 'rate of reaction');
model.param.set('PPD_MW_S', '108.144 [g/mol]', 'molar mass of substrate');
model.param.set('PPD_MW_P1', '150.181 [g/mol]', 'molar mass of product 1');
model.param.set('PPD_MW_P2', '192.218 [g/mol]', 'molar mass of product 2');
model.param.set('PPD_Kow_S', '10^-0.68', 'log partition coefficient of S in octanol/water');
model.param.set('PPD_Kow_P1', '10^0.08', 'log partition coefficient of P1 in octanol/water');
model.param.set('PPD_Kow_P2', '10^-0.45', 'log partition coefficient of P2 in octanol/water');
model.param.set('PPD_fu_S', '0.6625', 'fraction of unbound solute S');
model.param.set('PPD_fu_P1', '0.8102', 'fraction of unbound solute P1');
model.param.set('PPD_fu_P2', '0.7925', 'fraction of unbound solute P2');
model.param.set('PPD_fnon_S', '0.9264', 'fraction of un-ionised solute S');
model.param.set('PPD_fnon_P1', '0.998', 'fraction of un-ionised solute P1');
model.param.set('PPD_fnon_P2', '0.1', 'fraction of un-ionised solute P2');

% Diffusion coefficients using Longjian Chen's model
model.param.set('AHT_D_sc_S_fine', '1.33384*10^-13', 'AHT S diffusion coefficient, fine');
model.param.set('AHT_D_sc_S_coarse', '1.33289*10^-13', 'AHT S diffusion coefficient, coarse');
model.param.set('AHT_D_sc_P1_fine', '9.18883*10^-14', 'AHT P1 diffusion coefficient, fine');
model.param.set('AHT_D_sc_P1_coarse', '9.18203*10^-14', 'AHT P1 diffusion coefficient, coarse');
model.param.set('AHT_D_sc_P2_fine', '6.86547*10^-14', 'AHT P2 diffusion coefficient, fine');
model.param.set('AHT_D_sc_P2_coarse', '6.85890*10^-14', 'AHT P2 diffusion coefficient, coarse');
model.param.set('AHT_D_sc_P3_fine', '4.19145*10^-14', 'AHT P3 diffusion coefficient, fine');
model.param.set('AHT_D_sc_P3_coarse', '4.18822*10^-14', 'AHT P3 diffusion coefficient, coarse');
model.param.set('AMC_D_sc_S_fine', '1.27054*10^-13', 'AMC S diffusion coefficient, fine');
model.param.set('AMC_D_sc_S_coarse', '1.26944*10^-13', 'AMC S diffusion coefficient, coarse');
model.param.set('AMC_D_sc_P1_fine', '8.98866*10^-14', 'AMC P1 diffusion coefficient, fine');
model.param.set('AMC_D_sc_P1_coarse', '8.98166*10^-14', 'AMC P1 diffusion coefficient, coarse');
model.param.set('AMC_D_sc_P2_fine', '6.89537*10^-14', 'AMC P2 diffusion coefficient, fine');
model.param.set('AMC_D_sc_P2_coarse', '6.88922*10^-14', 'AMC P2 diffusion coefficient, coarse');
model.param.set('AMC_D_sc_P3_fine', '4.20832*10^-14', 'AMC P3 diffusion coefficient, fine');
model.param.set('AMC_D_sc_P3_coarse', '4.20513*10^-14', 'AMC P3 diffusion coefficient, coarse');
model.param.set('AEP_D_sc_S_fine', '1.23929*10^-13', 'AEP S diffusion coefficient, fine');
model.param.set('AEP_D_sc_S_coarse', '1.23850*10^-13', 'AEP S diffusion coefficient, coarse');
model.param.set('AEP_D_sc_P1_fine', '8.47459*10^-14', 'AEP P1 diffusion coefficient, fine');
model.param.set('AEP_D_sc_P1_coarse', '8.4688*10^-14', 'AEP P1 diffusion coefficient, coarse');
model.param.set('AEP_D_sc_P2_fine', '6.28108*10^-14', 'AEP P2 diffusion coefficient, fine');
model.param.set('AEP_D_sc_P2_coarse', '6.27553*10^-14', 'AEP P2 diffusion coefficient, coarse');
model.param.set('AEP_D_sc_P3_fine', '3.82388*10^-14', 'AEP P3 diffusion coefficient, fine');
model.param.set('AEP_D_sc_P3_coarse', '3.82040*10^-14', 'AEP P3 diffusion coefficient, coarse');
model.param.set('TDA_D_sc_S_fine', '1.26566*10^-13', 'TDA S diffusion coefficient, fine');
model.param.set('TDA_D_sc_S_coarse', '1.26453*10^-13', 'TDA S diffusion coefficient, coarse');
model.param.set('TDA_D_sc_P1_fine', '8.97779*10^-14', 'TDA P1 diffusion coefficient, fine');
model.param.set('TDA_D_sc_P1_coarse', '8.97033*10^-14', 'TDA P1 diffusion coefficient, coarse');
model.param.set('TDA_D_sc_P2_fine', '8.97779*10^-14', 'TDA P2 diffusion coefficient, fine');
model.param.set('TDA_D_sc_P2_coarse', '8.97033*10^-14', 'TDA P2 diffusion coefficient, coarse');
model.param.set('TDA_D_sc_P3_fine', '6.80119*10^-14', 'TDA P3 diffusion coefficient, fine');
model.param.set('TDA_D_sc_P3_coarse', '6.79485*10^-14', 'TDA P3 diffusion coefficient, coarse');
model.param.set('HDAP_D_sc_S_fine', '1.07725*10^-13', 'HDAP S diffusion coefficient, fine');
model.param.set('HDAP_D_sc_S_coarse', '1.07633*10^-13', 'HDAP S diffusion coefficient, coarse');
model.param.set('HDAP_D_sc_P1_fine', '8.07218*10^-14', 'HDAP P1 diffusion coefficient, fine');
model.param.set('HDAP_D_sc_P1_coarse', '8.06633*10^-14', 'HDAP P1 diffusion coefficient, coarse');
model.param.set('PPD_D_sc_S_fine', '1.44799*10^-13', 'PPD S diffusion coefficient, fine');
model.param.set('PPD_D_sc_S_coarse', '1.44664*10^-13', 'PPD S diffusion coefficient, coarse');
model.param.set('PPD_D_sc_P1_fine', '9.9388*10^-14', 'PPD P1 diffusion coefficient, fine');
model.param.set('PPD_D_sc_P1_coarse', '9.92957*10^-14', 'PPD P1 diffusion coefficient, coarse');
model.param.set('PPD_D_sc_P2_fine', '7.47247*10^-14', 'PPD P2 diffusion coefficient, fine');
model.param.set('PPD_D_sc_P2_coarse', '4.46547*10^-14', 'PPD P2 diffusion coefficient, coarse');

%% Unknown parameters for analysis
global vehic 
K_v_S = vehic %1; % Partition coefficient of PT, vehicle/water
global kr
alt_AHT_k_rate = kr*3.8055E-05
model.param.set('K_v_S', K_v_S, 'Partition coefficient of PT, vehicle/water');
model.param.set('k_rate', [num2str(alt_AHT_k_rate) '[1/s]'], 'rate of reaction, 3.8055E-05');

%% Components
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 1);

%% Geometry
model.component('comp1').geom('geom1').lengthUnit([native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').label('Vehicle');
model.component('comp1').geom('geom1').feature('i1').set('p2', 100);
model.component('comp1').geom('geom1').create('i2', 'Interval');
model.component('comp1').geom('geom1').feature('i2').label('Stratum corneum');
model.component('comp1').geom('geom1').feature('i2').set('p1', 100);
model.component('comp1').geom('geom1').feature('i2').set('p2', 110.5);
model.component('comp1').geom('geom1').create('i4', 'Interval');
model.component('comp1').geom('geom1').feature('i4').label('Viable epidermis');
model.component('comp1').geom('geom1').feature('i4').set('p1', 110.5);
model.component('comp1').geom('geom1').feature('i4').set('p2', 150.5);
model.component('comp1').geom('geom1').create('i5', 'Interval');
model.component('comp1').geom('geom1').feature('i5').label('Dermis');
model.component('comp1').geom('geom1').feature('i5').set('p1', 150.5);
model.component('comp1').geom('geom1').feature('i5').set('p2', 600);
model.component('comp1').geom('geom1').create('i3', 'Interval');
model.component('comp1').geom('geom1').feature('i3').label('Receptor');
model.component('comp1').geom('geom1').feature('i3').set('p1', 600);
model.component('comp1').geom('geom1').feature('i3').set('p2', 2100);
model.component('comp1').geom('geom1').run;

%% Definitions
% Integration
model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl.create('intop2', 'Integration');
model.component('comp1').cpl.create('intop3', 'Integration');
model.component('comp1').cpl.create('intop4', 'Integration');
model.component('comp1').cpl.create('intop5', 'Integration');
model.component('comp1').cpl('intop1').selection.set([1]);
model.component('comp1').cpl('intop2').selection.set([2]);
model.component('comp1').cpl('intop3').selection.set([3]);
model.component('comp1').cpl('intop4').selection.set([4]);
model.component('comp1').cpl('intop5').selection.set([5]);

%% Transport of dilute species
% Define Transport of dilute species - vehicle
model.component('comp1').physics.create('tds', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds').field('concentration').field('c1');
model.component('comp1').physics('tds').field('concentration').component({'c1'});
model.component('comp1').physics('tds').selection.set([1]);
model.component('comp1').physics('tds').create('constr1', 'PointwiseConstraint', 0);
model.component('comp1').physics('tds').feature('constr1').selection.set([2]);

% Define Transport of dilute species - stratum corneum
model.component('comp1').physics.create('tds2', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds2').field('concentration').field('c2');
model.component('comp1').physics('tds2').field('concentration').component({'c2'});
model.component('comp1').physics('tds2').selection.set([2]);
model.component('comp1').physics('tds2').create('constr1', 'PointwiseConstraint', 0);
model.component('comp1').physics('tds2').feature('constr1').selection.set([3]);

% Define Transport of dilute species - viable epidermis
model.component('comp1').physics.create('tds3', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds3').field('concentration').field('c3s');
model.component('comp1').physics('tds3').field('concentration').component({'c3s' 'c3p'});
model.component('comp1').physics('tds3').selection.set([3]);
model.component('comp1').physics('tds3').create('constr1', 'PointwiseConstraint', 0);
model.component('comp1').physics('tds3').feature('constr1').selection.set([4]);
model.component('comp1').physics('tds3').create('constr2', 'PointwiseConstraint', 0);
model.component('comp1').physics('tds3').feature('constr2').selection.set([4]);
model.component('comp1').physics('tds3').create('reac1', 'Reactions', 1);
model.component('comp1').physics('tds3').feature('reac1').selection.set([3]);

% Define Transport of dilute species - dermis
model.component('comp1').physics.create('tds4', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds4').field('concentration').field('c4s');
model.component('comp1').physics('tds4').field('concentration').component({'c4s' 'c4p'});
model.component('comp1').physics('tds4').selection.set([4]);
model.component('comp1').physics('tds4').create('constr1', 'PointwiseConstraint', 0);
model.component('comp1').physics('tds4').feature('constr1').selection.set([5]);
model.component('comp1').physics('tds4').create('constr2', 'PointwiseConstraint', 0);
model.component('comp1').physics('tds4').feature('constr2').selection.set([5]);

% Define Transport of dilute species - receptor fluid
model.component('comp1').physics.create('tds5', 'DilutedSpecies', 'geom1');
model.component('comp1').physics('tds5').field('concentration').field('c5s');
model.component('comp1').physics('tds5').field('concentration').component({'c5s' 'c5p'});
model.component('comp1').physics('tds5').selection.set([5]);

% Details of Transport of dilute species - vehicle
model.component('comp1').physics('tds').label('Vehicle');
model.component('comp1').physics('tds').prop('TransportMechanism').set('Convection', false);
model.component('comp1').physics('tds').feature('cdm1').set('minput_temperature', 'T');
model.component('comp1').physics('tds').feature('cdm1').set('D_c1', {'D_v_S'; '0'; '0'; '0'; 'D_v_S'; '0'; '0'; '0'; 'D_v_S'});
model.component('comp1').physics('tds').feature('init1').set('initc', 'conc');
model.component('comp1').physics('tds').feature('constr1').set('constraintType', 'userDefined');
model.component('comp1').physics('tds').feature('constr1').set('constraintExpression', '(t<3600)*(c1-(K_vsc_S*c2))');
model.component('comp1').physics('tds').feature('constr1').set('constraintForce', 'test(c1)-test(c2)');

% Details of Transport of dilute species - stratum corneum
model.component('comp1').physics('tds2').label('Stratum corneum');
model.component('comp1').physics('tds2').prop('TransportMechanism').set('Convection', false);
model.component('comp1').physics('tds2').feature('cdm1').set('D_c2', {'D_sc_S'; '0'; '0'; '0'; 'D_sc_S'; '0'; '0'; '0'; 'D_sc_S'});
model.component('comp1').physics('tds2').feature('cdm1').set('minput_temperature', 'T');
model.component('comp1').physics('tds2').feature('constr1').set('constraintType', 'userDefined');
model.component('comp1').physics('tds2').feature('constr1').set('constraintExpression', 'c2-(K_scve_S*c3s)');
model.component('comp1').physics('tds2').feature('constr1').set('constraintForce', 'test(c2)-test(c3s)');

% Details of Transport of dilute species - viable epidermis
model.component('comp1').physics('tds3').label('Viable epidermis');
model.component('comp1').physics('tds3').prop('TransportMechanism').set('Convection', false);
model.component('comp1').physics('tds3').feature('cdm1').set('minput_temperature', 'T');
model.component('comp1').physics('tds3').feature('cdm1').set('D_c3s', {'D_ve_S'; '0'; '0'; '0'; 'D_ve_S'; '0'; '0'; '0'; 'D_ve_S'});
model.component('comp1').physics('tds3').feature('cdm1').set('D_c3p', {'D_ve_P1'; '0'; '0'; '0'; 'D_ve_P1'; '0'; '0'; '0'; 'D_ve_P1'});
model.component('comp1').physics('tds3').feature('constr1').set('constraintType', 'userDefined');
model.component('comp1').physics('tds3').feature('constr1').set('constraintExpression', 'c3s-(K_ved_S*c4s)');
model.component('comp1').physics('tds3').feature('constr1').set('constraintForce', 'test(c3s)-test(c4s)');
model.component('comp1').physics('tds3').feature('constr2').set('constraintType', 'userDefined');
model.component('comp1').physics('tds3').feature('constr2').set('constraintExpression', 'c3p-(K_ved_P1*c4p)');
model.component('comp1').physics('tds3').feature('constr2').set('constraintForce', 'test(c3p)-test(c4p)');
model.component('comp1').physics('tds3').feature('reac1').set('R_c3s', '-k_rate*c3s');
model.component('comp1').physics('tds3').feature('reac1').set('R_c3p', 'k_rate*c3s');

% Details of Transport of dilute species - dermis
model.component('comp1').physics('tds4').label('Dermis');
model.component('comp1').physics('tds4').prop('TransportMechanism').set('Convection', false);
model.component('comp1').physics('tds4').feature('cdm1').set('minput_temperature', 'T');
model.component('comp1').physics('tds4').feature('cdm1').set('D_c4s', {'D_d_S'; '0'; '0'; '0'; 'D_d_S'; '0'; '0'; '0'; 'D_d_S'});
model.component('comp1').physics('tds4').feature('cdm1').set('D_c4p', {'D_d_P1'; '0'; '0'; '0'; 'D_d_P1'; '0'; '0'; '0'; 'D_d_P1'});
model.component('comp1').physics('tds4').feature('constr1').set('constraintType', 'userDefined');
model.component('comp1').physics('tds4').feature('constr1').set('constraintExpression', 'c4s-(K_dr_S*c5s)');
model.component('comp1').physics('tds4').feature('constr1').set('constraintForce', 'test(c4s)-test(c5s)');
model.component('comp1').physics('tds4').feature('constr2').set('constraintType', 'userDefined');
model.component('comp1').physics('tds4').feature('constr2').set('constraintExpression', 'c4p-(K_dr_P1*c5p)');
model.component('comp1').physics('tds4').feature('constr2').set('constraintForce', 'test(c4p)-test(c5p)');

% Details of Transport of dilute species - receptor fluid
model.component('comp1').physics('tds5').label('Receptor');
model.component('comp1').physics('tds5').prop('TransportMechanism').set('Convection', false);
model.component('comp1').physics('tds5').prop('MassConsistentStabilization').set('glim_mass', '(0.1[mol/m^3])/tds3.helem');
model.component('comp1').physics('tds5').feature('cdm1').set('minput_temperature', 'T');
model.component('comp1').physics('tds5').feature('cdm1').set('D_c5s', {'D_r_S'; '0'; '0'; '0'; 'D_r_S'; '0'; '0'; '0'; 'D_r_S'});
model.component('comp1').physics('tds5').feature('cdm1').set('D_c5p', {'D_r_P1'; '0'; '0'; '0'; 'D_r_P1'; '0'; '0'; '0'; 'D_r_P1'});

%% Mesh
% Create mesh
model.component('comp1').mesh.create('mesh1');

% Mesh selection
model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').create('size2', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').create('size3', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').selection.geom('geom1', 1);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').selection.geom('geom1', 1);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').selection.set([1 3 4]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size3').selection.geom('geom1', 1);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size3').selection.set([5]);

% Mesh definition
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', 5);
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 0.102);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').label('Size 1 - SC');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 0.05);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmin', 1.53);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').label('Size 2 - Everything else');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hmax', 1);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hmin', 1.53);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size3').label('Size 2 - Everything else 1');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size3').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size3').set('hmax', 500);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size3').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size3').set('hmin', 1.53);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size3').set('hminactive', false);
model.component('comp1').mesh('mesh1').run;
model.frame('mesh1').sorder(1);

%% Study
% Create study
model.study.create('std1');
model.study('std1').create('time', 'Transient');
% Create solution
model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature.remove('fcDef');
%
model.result.numerical.create('int1', 'IntLine');
model.result.numerical.create('int2', 'IntLine');
model.result.numerical.create('int3', 'IntLine');
model.result.numerical('int1').set('probetag', 'none');
model.result.numerical('int2').set('probetag', 'none');
model.result.numerical('int3').set('probetag', 'none');

%% Create plots
% Concentration plots
model.result.create('pg1', 'PlotGroup1D');
model.result.create('pg2', 'PlotGroup1D');
model.result.create('pg3', 'PlotGroup1D');
model.result.create('pg8', 'PlotGroup1D');
model.result.create('pg9', 'PlotGroup1D');
model.result.create('pg10', 'PlotGroup1D');
model.result.create('pg11', 'PlotGroup1D');
model.result.create('pg12', 'PlotGroup1D');
% Integration plots
model.result.create('pg4', 'PlotGroup1D');
model.result.create('pg5', 'PlotGroup1D');
model.result.create('pg6', 'PlotGroup1D');
model.result.create('pg13', 'PlotGroup1D');
model.result.create('pg14', 'PlotGroup1D');
model.result.create('pg15', 'PlotGroup1D');
model.result.create('pg16', 'PlotGroup1D');
model.result.create('pg17', 'PlotGroup1D');
% Total integration plot (mass balance)
model.result.create('pg7', 'PlotGroup1D');
% Concentration plots hourly
model.result.create('pg18', 'PlotGroup1D');
model.result.create('pg19', 'PlotGroup1D');
model.result.create('pg20', 'PlotGroup1D');
model.result.create('pg21', 'PlotGroup1D');
model.result.create('pg22', 'PlotGroup1D');
model.result.create('pg23', 'PlotGroup1D');
model.result.create('pg25', 'PlotGroup1D');
model.result.create('pg24', 'PlotGroup1D');
% Integration plots end point
model.result.create('pg26', 'PlotGroup1D');
model.result.create('pg27', 'PlotGroup1D');
model.result.create('pg28', 'PlotGroup1D');
model.result.create('pg29', 'PlotGroup1D');
model.result.create('pg30', 'PlotGroup1D');
model.result.create('pg31', 'PlotGroup1D');
model.result.create('pg32', 'PlotGroup1D');
model.result.create('pg33', 'PlotGroup1D');

% Define data range to be plotted
% Concentration plots
model.result('pg1').create('lngr1', 'LineGraph');
model.result('pg1').feature('lngr1').selection.set([1]);
model.result('pg2').create('lngr1', 'LineGraph');
model.result('pg2').feature('lngr1').selection.set([2]);
model.result('pg3').create('lngr1', 'LineGraph');
model.result('pg3').feature('lngr1').selection.set([3]);
model.result('pg8').create('lngr1', 'LineGraph');
model.result('pg8').feature('lngr1').selection.set([3]);
model.result('pg9').create('lngr1', 'LineGraph');
model.result('pg9').feature('lngr1').selection.set([4]);
model.result('pg10').create('lngr1', 'LineGraph');
model.result('pg10').feature('lngr1').selection.set([4]);
model.result('pg11').create('lngr1', 'LineGraph');
model.result('pg11').feature('lngr1').selection.set([5]);
model.result('pg12').create('lngr1', 'LineGraph');
model.result('pg12').feature('lngr1').selection.set([5]);
% Integration plots
model.result('pg4').create('ptgr1', 'PointGraph');
model.result('pg4').feature('ptgr1').set('data', 'dset1');
model.result('pg4').feature('ptgr1').selection.set([1 2]);
model.result('pg5').create('ptgr1', 'PointGraph');
model.result('pg5').feature('ptgr1').set('data', 'dset1');
model.result('pg5').feature('ptgr1').selection.set([2 3]);
model.result('pg6').create('ptgr1', 'PointGraph');
model.result('pg6').feature('ptgr1').set('data', 'dset1');
model.result('pg6').feature('ptgr1').selection.set([3 4]);
model.result('pg13').create('ptgr1', 'PointGraph');
model.result('pg13').feature('ptgr1').set('data', 'dset1');
model.result('pg13').feature('ptgr1').selection.set([3 4]);
model.result('pg14').create('ptgr1', 'PointGraph');
model.result('pg14').feature('ptgr1').set('data', 'dset1');
model.result('pg14').feature('ptgr1').selection.set([4 5]);
model.result('pg15').create('ptgr1', 'PointGraph');
model.result('pg15').feature('ptgr1').set('data', 'dset1');
model.result('pg15').feature('ptgr1').selection.set([4 5]);
model.result('pg16').create('ptgr1', 'PointGraph');
model.result('pg16').feature('ptgr1').set('data', 'dset1');
model.result('pg16').feature('ptgr1').selection.set([5 6]);
model.result('pg17').create('ptgr1', 'PointGraph');
model.result('pg17').feature('ptgr1').set('data', 'dset1');
model.result('pg17').feature('ptgr1').selection.set([5 6]);
% Total integration plot (mass balance)
model.result('pg7').create('ptgr1', 'PointGraph');
model.result('pg7').feature('ptgr1').selection.all;
% Concentration plots hourly
model.result('pg18').create('lngr1', 'LineGraph');
model.result('pg18').feature('lngr1').selection.set([1]);
model.result('pg19').create('lngr1', 'LineGraph');
model.result('pg19').feature('lngr1').selection.set([2]);
model.result('pg20').create('lngr1', 'LineGraph');
model.result('pg20').feature('lngr1').selection.set([3]);
model.result('pg21').create('lngr1', 'LineGraph');
model.result('pg21').feature('lngr1').selection.set([3]);
model.result('pg22').create('lngr1', 'LineGraph');
model.result('pg22').feature('lngr1').selection.set([4]);
model.result('pg23').create('lngr1', 'LineGraph');
model.result('pg23').feature('lngr1').selection.set([4]);
model.result('pg25').create('lngr1', 'LineGraph');
model.result('pg25').feature('lngr1').selection.set([5]);
model.result('pg24').create('lngr1', 'LineGraph');
model.result('pg24').feature('lngr1').selection.set([5]);
% Integration plots end points
model.result('pg26').create('ptgr1', 'PointGraph');
model.result('pg26').feature('ptgr1').set('data', 'dset1');
model.result('pg26').feature('ptgr1').selection.set([1 2]);
model.result('pg27').create('ptgr1', 'PointGraph');
model.result('pg27').feature('ptgr1').set('data', 'dset1');
model.result('pg27').feature('ptgr1').selection.set([2 3]);
model.result('pg28').create('ptgr1', 'PointGraph');
model.result('pg28').feature('ptgr1').set('data', 'dset1');
model.result('pg28').feature('ptgr1').selection.set([3 4]);
model.result('pg29').create('ptgr1', 'PointGraph');
model.result('pg29').feature('ptgr1').set('data', 'dset1');
model.result('pg29').feature('ptgr1').selection.set([3 4]);
model.result('pg30').create('ptgr1', 'PointGraph');
model.result('pg30').feature('ptgr1').set('data', 'dset1');
model.result('pg30').feature('ptgr1').selection.set([4 5]);
model.result('pg31').create('ptgr1', 'PointGraph');
model.result('pg31').feature('ptgr1').set('data', 'dset1');
model.result('pg31').feature('ptgr1').selection.set([4 5]);
model.result('pg32').create('ptgr1', 'PointGraph');
model.result('pg32').feature('ptgr1').set('data', 'dset1');
model.result('pg32').feature('ptgr1').selection.set([5 6]);
model.result('pg33').create('ptgr1', 'PointGraph');
model.result('pg33').feature('ptgr1').set('data', 'dset1');
model.result('pg33').feature('ptgr1').selection.set([5 6]);

%% Transient
model.study('std1').feature('time').set('tlist', 'range(0,1,86400)');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('v1').set('clist', {'range(0,1,86400)'});
model.sol('sol1').feature('t1').set('tlist', 'range(0,1,86400)');
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 8);
model.sol('sol1').feature('t1').feature('fc1').set('damp', 0.9);
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').runAll;

%% Label result plots
%Concentration plots
model.result('pg1').label('Concentration (tds)');
model.result('pg1').set('xlabel', 'Arc length');
model.result('pg1').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg1').set('ylog', true);
model.result('pg1').set('xlabelactive', false);
model.result('pg1').set('ylabelactive', false);
model.result('pg1').feature('lngr1').label('Line Graph');
model.result('pg1').feature('lngr1').set('resolution', 'normal');
model.result('pg2').label('Concentration (tds2)');
model.result('pg2').set('xlabel', 'Arc length');
model.result('pg2').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg2').set('xlabelactive', false);
model.result('pg2').set('ylabelactive', false);
model.result('pg2').feature('lngr1').label('Line Graph');
model.result('pg2').feature('lngr1').set('expr', 'c2');
model.result('pg2').feature('lngr1').set('resolution', 'normal');
model.result('pg3').label('Concentration (tds3)');
model.result('pg3').set('xlabel', 'Arc length');
model.result('pg3').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg3').set('xlabelactive', false);
model.result('pg3').set('ylabelactive', false);
model.result('pg3').feature('lngr1').label('Line Graph');
model.result('pg3').feature('lngr1').set('expr', 'c3s');
model.result('pg3').feature('lngr1').set('resolution', 'normal');
model.result('pg8').label('Concentration (tds4) ');
model.result('pg8').set('xlabel', 'Arc length');
model.result('pg8').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg8').set('xlabelactive', false);
model.result('pg8').set('ylabelactive', false);
model.result('pg8').feature('lngr1').set('expr', 'c3p');
model.result('pg8').feature('lngr1').set('resolution', 'normal');
model.result('pg9').label('Concentration (tds5)');
model.result('pg9').set('xlabel', 'Arc length');
model.result('pg9').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg9').set('xlabelactive', false);
model.result('pg9').set('ylabelactive', false);
model.result('pg9').feature('lngr1').set('expr', 'c4s');
model.result('pg9').feature('lngr1').set('resolution', 'normal');
model.result('pg10').label('Concentration (tds6)');
model.result('pg10').set('xlabel', 'Arc length');
model.result('pg10').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg10').set('xlabelactive', false);
model.result('pg10').set('ylabelactive', false);
model.result('pg10').feature('lngr1').set('expr', 'c4p');
model.result('pg10').feature('lngr1').set('resolution', 'normal');
model.result('pg11').label('Concentration (tds7)');
model.result('pg11').set('xlabel', 'Arc length');
model.result('pg11').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg11').set('xlabelactive', false);
model.result('pg11').set('ylabelactive', false);
model.result('pg11').feature('lngr1').set('expr', 'c5s');
model.result('pg11').feature('lngr1').set('resolution', 'normal');
model.result('pg12').label('Concentration (tds8)');
model.result('pg12').set('xlabel', 'Arc length');
model.result('pg12').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg12').set('xlabelactive', false);
model.result('pg12').set('ylabelactive', false);
model.result('pg12').feature('lngr1').set('expr', 'c5p');
model.result('pg12').feature('lngr1').set('resolution', 'normal');
% Integration plots
model.result('pg4').set('xlabel', 'Time (s)');
model.result('pg4').set('ylabel', 'Integration 1 (mol/m<sup>2</sup>)');
model.result('pg4').set('xlabelactive', false);
model.result('pg4').set('ylabelactive', false);
model.result('pg4').feature('ptgr1').set('expr', 'intop1(c1)');
model.result('pg4').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg4').feature('ptgr1').set('descr', 'Integration 1');
model.result('pg5').set('xlabel', 'Time (s)');
model.result('pg5').set('ylabel', 'Integration 2 (mol/m<sup>2</sup>)');
model.result('pg5').set('xlabelactive', false);
model.result('pg5').set('ylabelactive', false);
model.result('pg5').feature('ptgr1').set('expr', 'intop2(c2)');
model.result('pg5').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg5').feature('ptgr1').set('descr', 'Integration 2');
model.result('pg6').set('xlabel', 'Time (s)');
model.result('pg6').set('ylabel', 'Integration 5 (mol/m<sup>2</sup>)');
model.result('pg6').set('xlabelactive', false);
model.result('pg6').set('ylabelactive', false);
model.result('pg6').feature('ptgr1').set('expr', 'intop3(c3s)');
model.result('pg6').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg6').feature('ptgr1').set('descr', 'Integration 3');
model.result('pg13').set('xlabel', 'Time (s)');
model.result('pg13').set('ylabel', 'Integration 3 (mol/m<sup>2</sup>)');
model.result('pg13').set('xlabelactive', false);
model.result('pg13').set('ylabelactive', false);
model.result('pg13').feature('ptgr1').set('expr', 'intop3(c3p)');
model.result('pg13').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg13').feature('ptgr1').set('descr', 'Integration 3');
model.result('pg14').set('xlabel', 'Time (s)');
model.result('pg14').set('ylabel', 'Integration 4 (mol/m<sup>2</sup>)');
model.result('pg14').set('xlabelactive', false);
model.result('pg14').set('ylabelactive', false);
model.result('pg14').feature('ptgr1').set('expr', 'intop4(c4s)');
model.result('pg14').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg14').feature('ptgr1').set('descr', 'Integration 4');
model.result('pg15').set('xlabel', 'Time (s)');
model.result('pg15').set('ylabel', 'Integration 4 (mol/m<sup>2</sup>)');
model.result('pg15').set('xlabelactive', false);
model.result('pg15').set('ylabelactive', false);
model.result('pg15').feature('ptgr1').set('expr', 'intop4(c4p)');
model.result('pg15').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg15').feature('ptgr1').set('descr', 'Integration 4');
model.result('pg16').set('xlabel', 'Time (s)');
model.result('pg16').set('ylabel', 'Integration 5 (mol/m<sup>2</sup>)');
model.result('pg16').set('xlabelactive', false);
model.result('pg16').set('ylabelactive', false);
model.result('pg16').feature('ptgr1').set('expr', 'intop5(c5s)');
model.result('pg16').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg16').feature('ptgr1').set('descr', 'Integration 5');
model.result('pg17').set('xlabel', 'Time (s)');
model.result('pg17').set('ylabel', 'Integration 5 (mol/m<sup>2</sup>)');
model.result('pg17').set('xlabelactive', false);
model.result('pg17').set('ylabelactive', false);
model.result('pg17').feature('ptgr1').set('expr', 'intop5(c5p)');
model.result('pg17').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg17').feature('ptgr1').set('descr', 'Integration 5');
% Total integration plot (mass balance)
model.result('pg7').set('xlabel', 'Time (s)');
model.result('pg7').set('ylabel', 'intop1(c1)+intop2(c2)+intop3(c3s)+intop3(c3p)+intop4(c4s)+intop4(c4p)+intop5(c5s)+intop5(c5p) (mol/m<sup>2</sup>)');
model.result('pg7').set('xlabelactive', false);
model.result('pg7').set('ylabelactive', false);
model.result('pg7').feature('ptgr1').set('expr', 'intop1(c1)+intop2(c2)+intop3(c3s)+intop3(c3p)+intop4(c4s)+intop4(c4p)+intop5(c5s)+intop5(c5p)');
model.result('pg7').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg7').feature('ptgr1').set('descr', 'intop1(c1)+intop2(c2)+intop3(c3s)+intop3(c3p)+intop4(c4s)+intop4(c4p)+intop5(c5s)+intop5(c5p)');
% Concentration plots hourly
model.result('pg18').set('looplevelinput', {'manualindices'});
model.result('pg18').set('looplevelindices', {'range(1,3600,86401)'});
model.result('pg18').set('xlabel', 'Arc length');
model.result('pg18').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg18').set('xlabelactive', false);
model.result('pg18').set('ylabelactive', false);
model.result('pg18').feature('lngr1').set('resolution', 'normal');
model.result('pg19').set('looplevelinput', {'manualindices'});
model.result('pg19').set('looplevelindices', {'range(1,3600,86401)'});
model.result('pg19').set('xlabel', 'Arc length');
model.result('pg19').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg19').set('xlabelactive', false);
model.result('pg19').set('ylabelactive', false);
model.result('pg19').feature('lngr1').set('expr', 'c2');
model.result('pg19').feature('lngr1').set('resolution', 'normal');
model.result('pg20').set('looplevelinput', {'manualindices'});
model.result('pg20').set('looplevelindices', {'range(1,3600,86401)'});
model.result('pg20').set('xlabel', 'Arc length');
model.result('pg20').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg20').set('xlabelactive', false);
model.result('pg20').set('ylabelactive', false);
model.result('pg20').feature('lngr1').set('expr', 'c3s');
model.result('pg20').feature('lngr1').set('resolution', 'normal');
model.result('pg21').set('looplevelinput', {'manualindices'});
model.result('pg21').set('looplevelindices', {'range(1,3600,86401)'});
model.result('pg21').set('xlabel', 'Arc length');
model.result('pg21').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg21').set('xlabelactive', false);
model.result('pg21').set('ylabelactive', false);
model.result('pg21').feature('lngr1').set('expr', 'c3p');
model.result('pg21').feature('lngr1').set('resolution', 'normal');
model.result('pg22').set('looplevelinput', {'manualindices'});
model.result('pg22').set('looplevelindices', {'range(1,3600,86401)'});
model.result('pg22').set('xlabel', 'Arc length');
model.result('pg22').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg22').set('xlabelactive', false);
model.result('pg22').set('ylabelactive', false);
model.result('pg22').feature('lngr1').set('expr', 'c4s');
model.result('pg22').feature('lngr1').set('resolution', 'normal');
model.result('pg23').set('looplevelinput', {'manualindices'});
model.result('pg23').set('looplevelindices', {'range(1,3600,86401)'});
model.result('pg23').set('xlabel', 'Arc length');
model.result('pg23').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg23').set('xlabelactive', false);
model.result('pg23').set('ylabelactive', false);
model.result('pg23').feature('lngr1').set('expr', 'c4p');
model.result('pg23').feature('lngr1').set('resolution', 'normal');
model.result('pg25').set('looplevelinput', {'manualindices'});
model.result('pg25').set('looplevelindices', {'range(1,3600,86401)'});
model.result('pg25').set('xlabel', 'Arc length');
model.result('pg25').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg25').set('xlabelactive', false);
model.result('pg25').set('ylabelactive', false);
model.result('pg25').feature('lngr1').set('expr', 'c5s');
model.result('pg25').feature('lngr1').set('resolution', 'normal');
model.result('pg24').set('looplevelinput', {'manualindices'});
model.result('pg24').set('looplevelindices', {'range(1,3600,86401)'});
model.result('pg24').set('xlabel', 'Arc length');
model.result('pg24').set('ylabel', 'Concentration (mol/m<sup>3</sup>)');
model.result('pg24').set('xlabelactive', false);
model.result('pg24').set('ylabelactive', false);
model.result('pg24').feature('lngr1').set('expr', 'c5p');
model.result('pg24').feature('lngr1').set('resolution', 'normal');
% Integration plots end points
model.result('pg26').set('xlabel', 'Time (s)');
model.result('pg26').set('ylabel', 'Integration 1 (mol/m<sup>2</sup>)');
model.result('pg26').set('xlabelactive', false);
model.result('pg26').set('ylabelactive', false);
model.result('pg26').feature('ptgr1').set('looplevelinput', {'manualindices'});
model.result('pg26').feature('ptgr1').set('looplevelindices', [86401]);
model.result('pg26').feature('ptgr1').set('expr', 'intop1(c1)');
model.result('pg26').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg26').feature('ptgr1').set('descr', 'Integration 1');
model.result('pg27').set('xlabel', 'Time (s)');
model.result('pg27').set('ylabel', 'Integration 2 (mol/m<sup>2</sup>)');
model.result('pg27').set('xlabelactive', false);
model.result('pg27').set('ylabelactive', false);
model.result('pg27').feature('ptgr1').set('looplevelinput', {'manualindices'});
model.result('pg27').feature('ptgr1').set('looplevelindices', [86401]);
model.result('pg27').feature('ptgr1').set('expr', 'intop2(c2)');
model.result('pg27').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg27').feature('ptgr1').set('descr', 'Integration 2');
model.result('pg28').set('xlabel', 'Time (s)');
model.result('pg28').set('ylabel', 'Integration 5 (mol/m<sup>2</sup>)');
model.result('pg28').set('xlabelactive', false);
model.result('pg28').set('ylabelactive', false);
model.result('pg28').feature('ptgr1').set('looplevelinput', {'manualindices'});
model.result('pg28').feature('ptgr1').set('looplevelindices', [86401]);
model.result('pg28').feature('ptgr1').set('expr', 'intop3(c3s)');
model.result('pg28').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg28').feature('ptgr1').set('descr', 'Integration 3');
model.result('pg29').set('xlabel', 'Time (s)');
model.result('pg29').set('ylabel', 'Integration 3 (mol/m<sup>2</sup>)');
model.result('pg29').set('xlabelactive', false);
model.result('pg29').set('ylabelactive', false);
model.result('pg29').feature('ptgr1').set('looplevelinput', {'manualindices'});
model.result('pg29').feature('ptgr1').set('looplevelindices', [86401]);
model.result('pg29').feature('ptgr1').set('expr', 'intop3(c3p)');
model.result('pg29').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg29').feature('ptgr1').set('descr', 'Integration 3');
model.result('pg30').set('xlabel', 'Time (s)');
model.result('pg30').set('ylabel', 'Integration 4 (mol/m<sup>2</sup>)');
model.result('pg30').set('xlabelactive', false);
model.result('pg30').set('ylabelactive', false);
model.result('pg30').feature('ptgr1').set('looplevelinput', {'manualindices'});
model.result('pg30').feature('ptgr1').set('looplevelindices', [86401]);
model.result('pg30').feature('ptgr1').set('expr', 'intop4(c4s)');
model.result('pg30').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg30').feature('ptgr1').set('descr', 'Integration 4');
model.result('pg31').set('xlabel', 'Time (s)');
model.result('pg31').set('ylabel', 'Integration 4 (mol/m<sup>2</sup>)');
model.result('pg31').set('xlabelactive', false);
model.result('pg31').set('ylabelactive', false);
model.result('pg31').feature('ptgr1').set('looplevelinput', {'manualindices'});
model.result('pg31').feature('ptgr1').set('looplevelindices', [86401]);
model.result('pg31').feature('ptgr1').set('expr', 'intop4(c4p)');
model.result('pg31').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg31').feature('ptgr1').set('descr', 'Integration 4');
model.result('pg32').set('xlabel', 'Time (s)');
model.result('pg32').set('ylabel', 'Integration 5 (mol/m<sup>2</sup>)');
model.result('pg32').set('xlabelactive', false);
model.result('pg32').set('ylabelactive', false);
model.result('pg32').feature('ptgr1').set('looplevelinput', {'manualindices'});
model.result('pg32').feature('ptgr1').set('looplevelindices', [86401]);
model.result('pg32').feature('ptgr1').set('expr', 'intop5(c5s)');
model.result('pg32').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg32').feature('ptgr1').set('descr', 'Integration 5');
model.result('pg32').set('xlabel', 'Time (s)');
model.result('pg33').set('ylabel', 'Integration 5 (mol/m<sup>2</sup>)');
model.result('pg33').set('xlabelactive', false);
model.result('pg33').set('ylabelactive', false);
model.result('pg33').feature('ptgr1').set('looplevelinput', {'manualindices'});
model.result('pg33').feature('ptgr1').set('looplevelindices', [86401]);
model.result('pg33').feature('ptgr1').set('expr', 'intop5(c5p)');
model.result('pg33').feature('ptgr1').set('unit', 'mol/m^2');
model.result('pg33').feature('ptgr1').set('descr', 'Integration 5');

%% Export 
% Create export plots
% Integration plots
%model.result.export.create('plot1', 'Plot');
%model.result.export.create('plot2', 'Plot');
%model.result.export.create('plot3', 'Plot');
%model.result.export.create('plot5', 'Plot');
%model.result.export.create('plot6', 'Plot');
%model.result.export.create('plot7', 'Plot');
%model.result.export.create('plot8', 'Plot');
%model.result.export.create('plot9', 'Plot');
model.result.export.create('plot4', 'Plot');
% Integration plots end points
model.result.export.create('plot10', 'Plot');
model.result.export.create('plot11', 'Plot');
model.result.export.create('plot12', 'Plot');
model.result.export.create('plot13', 'Plot');
model.result.export.create('plot14', 'Plot');
model.result.export.create('plot15', 'Plot');
model.result.export.create('plot16', 'Plot');
model.result.export.create('plot17', 'Plot');
% Define data to be plotted and destination
% Integration plots
%model.result.export('plot1').set('plotgroup', 'pg4');
%model.result.export('plot1').set('plot', 'ptgr1');
%model.result.export('plot1').set('filename', ['D:\test\AHT_v_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
%model.result.export('plot2').set('plotgroup', 'pg5');
%model.result.export('plot2').set('plot', 'ptgr1');
%model.result.export('plot2').set('filename', ['D:\test\AHT_sc_S_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
%model.result.export('plot3').set('plotgroup', 'pg6');
%model.result.export('plot3').set('plot', 'ptgr1');
%model.result.export('plot3').set('filename', ['D:\test\AHT_ve_S_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
%model.result.export('plot5').set('plotgroup', 'pg13');
%model.result.export('plot5').set('plot', 'ptgr1');
%model.result.export('plot5').set('filename', ['D:\test\AHT_ve_P_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
%model.result.export('plot6').set('plotgroup', 'pg14');
%model.result.export('plot6').set('plot', 'ptgr1');
%model.result.export('plot6').set('filename', ['D:\test\AHT_d_S_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
%model.result.export('plot7').set('plotgroup', 'pg15');
%model.result.export('plot7').set('plot', 'ptgr1');
%model.result.export('plot7').set('filename', ['D:\test\AHT_d_P_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
%model.result.export('plot8').set('plotgroup', 'pg16');
%model.result.export('plot8').set('plot', 'ptgr1');
%model.result.export('plot8').set('filename', ['D:\test\AHT_r_S_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
%model.result.export('plot9').set('plotgroup', 'pg17');
%model.result.export('plot9').set('plot', 'ptgr1');
%model.result.export('plot9').set('filename', ['D:\test\AHT_r_P_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
model.result.export('plot4').set('plotgroup', 'pg7');
model.result.export('plot4').set('plot', 'ptgr1');
model.result.export('plot4').set('filename', ['D:\test\AHT_t_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
% Integration plots end points
model.result.export('plot10').set('plotgroup', 'pg26');
model.result.export('plot10').set('plot', 'ptgr1');
model.result.export('plot10').set('filename', ['D:\test\ep_AHT_v_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
model.result.export('plot11').set('plotgroup', 'pg27');
model.result.export('plot11').set('plot', 'ptgr1');
model.result.export('plot11').set('filename', ['D:\test\ep_AHT_sc_S_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
model.result.export('plot12').set('plotgroup', 'pg28');
model.result.export('plot12').set('plot', 'ptgr1');
model.result.export('plot12').set('filename', ['D:\test\ep_AHT_ve_S_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
model.result.export('plot13').set('plotgroup', 'pg29');
model.result.export('plot13').set('plot', 'ptgr1');
model.result.export('plot13').set('filename', ['D:\test\ep_AHT_ve_P_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
model.result.export('plot14').set('plotgroup', 'pg30');
model.result.export('plot14').set('plot', 'ptgr1');
model.result.export('plot14').set('filename', ['D:\test\ep_AHT_d_S_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
model.result.export('plot15').set('plotgroup', 'pg31');
model.result.export('plot15').set('plot', 'ptgr1');
model.result.export('plot15').set('filename', ['D:\test\ep_AHT_d_P_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
model.result.export('plot16').set('plotgroup', 'pg32');
model.result.export('plot16').set('plot', 'ptgr1');
model.result.export('plot16').set('filename', ['D:\test\ep_AHT_r_S_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);
model.result.export('plot17').set('plotgroup', 'pg33');
model.result.export('plot17').set('plot', 'ptgr1');
model.result.export('plot17').set('filename', ['D:\test\ep_AHT_r_P_Kv_' num2str(vehic) '_rf_' num2str(rf_length) '_krate_' num2str(kr) '_Kr_' num2str(recept) '.txt']);

% Run export
%model.result.export('plot1').run;
%model.result.export('plot2').run;
%model.result.export('plot3').run;
%model.result.export('plot5').run;
%model.result.export('plot6').run;
%model.result.export('plot7').run;
%model.result.export('plot8').run;
%model.result.export('plot9').run;
model.result.export('plot4').run;
model.result.export('plot10').run;
model.result.export('plot11').run;
model.result.export('plot12').run;
model.result.export('plot13').run;
model.result.export('plot14').run;
model.result.export('plot15').run;
model.result.export('plot16').run;
model.result.export('plot17').run;

out = model;
