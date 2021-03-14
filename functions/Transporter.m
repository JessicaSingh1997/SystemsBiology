function checkTransporter = Transporter('y')
x = find(strcmp(model.mets,'y'))
positions = find(model.S(x,:))
checkTransporter = [model.rxns(positions) getRxn_cobraFormat(model,positions) num2cell(model.lb(positions)) num2cell(model.ub(positions))]