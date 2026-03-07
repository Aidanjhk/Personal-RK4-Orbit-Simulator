function dragCoefficient = dragCoefficientCalculator(R, Tatm,  M, vSatMagnitude, sigma_t, sigma_n, betaAngle, pSat)
    Vw = sqrt((pi*R*Tatm)/(2*M));
    
    S = vSatMagnitude/Vw;
    if norm(pSat,2) < 120e3
        dragCoefficient = (2*sigma_n(betaAngle)*(Vw/vSatMagnitude)*(sin(betaAngle)^2) + (2/(sqrt(pi)*(S^2)))*((2-sigma_n(betaAngle))*(sin(betaAngle)^2) + sigma_t(betaAngle)*(cos(betaAngle)^2))*exp(-(S^2)*(sin(betaAngle)^2)) ...
        + 2*((2-sigma_n(betaAngle))*((sin(betaAngle)^2)+(1/(2*(S^2))))+sigma_t(betaAngle)*(cos(betaAngle)^2))*sin(betaAngle)*erf(S*sin(betaAngle)));
    else
        dragCoefficient = 2*(sigma_t(betaAngle) + sigma_n(betaAngle)*(Vw/vSatMagnitude)*sin(betaAngle) +  (2 - sigma_n(betaAngle) - sigma_t(betaAngle))*(sin(betaAngle)^2))*sin(betaAngle);
    end
end