function y = gaussianConvolution(xin, yin, x, dx)
    % Convolve the input function yin(xin) with a Gaussian resolution
    % function.
    %
    % For each output point xo = x(kout), the calculation is
    %
    %   y(xo) = integral[ yin(x) * G(xo-x; sigma) dx ]
    %
    % where
    %
    %   G(xo-x; sigma) = 1/(sqrt(2*pi)*sigma)
    %                    * exp(-(xo-x)^2/(2*sigma^2))
    %
    % The input function yin(xin) is assumed to be piecewise linear
    % between the supplied xin points. The convolution integral over
    % each linear segment can then be evaluated analytically.
    %
    % sigma = dx(kout) is the Gaussian resolution width associated
    % with each output point.

    y = zeros(size(x));
    Nin = numel(xin);
    Nout = numel(x);

    % log(0.001). This is used to truncate the Gaussian when its
    % amplitude has fallen to 0.1% of its maximum.
    LOGRESLIMIT = -6.90775527898213703123;

    % Index of the input point near the left-hand edge of the
    % Gaussian integration range. This is carried between output
    % points to avoid repeatedly searching from the beginning.
    kin = 1;

    for kout = 1:Nout

        % Gaussian resolution width for this output point.
        sigma = dx(kout);

        % Output coordinate at which we want the convolved value.
        xo = x(kout);

        % Truncate the Gaussian when it has fallen to 0.1% of
        % its peak value.
        %
        % exp(-limit^2/(2*sigma^2)) = 0.001
        %
        % Therefore
        %
        %   limit = sqrt(-2*sigma^2*log(0.001))
        %         ~= 3.717*sigma
        %
        % so the convolution only needs to consider approximately
        %
        %   xo - 3.717*sigma  <= x <=  xo + 3.717*sigma.
        limit = sqrt(-2.0 * sigma * sigma * LOGRESLIMIT);


        % Find the first input point at or just after the left-hand
        % edge of the Gaussian integration range.
        while kin < Nin && xin(kin) < xo - limit
            kin = kin + 1;
        end

        % Move back one point if necessary so that kin is the input
        % point immediately before (or near) xo-limit.
        while kin > 1 && xin(kin) > xo - limit
            kin = kin - 1;
        end


        if sigma > 0

            % Perform the actual Gaussian convolution at xo.
            %
            % convolveGaussianPoint treats yin(xin) as a piecewise
            % linear function and integrates each linear segment
            % analytically against the Gaussian resolution function.
            y(kout) = convolveGaussianPoint( ...
                xin, yin, kin, Nin, xo, limit, sigma);

        elseif kin < Nin

            % If sigma = 0 there is no resolution broadening.
            % Simply linearly interpolate yin at xo.

            m = (yin(kin + 1) - yin(kin)) / ...
                (xin(kin + 1) - xin(kin));

            b = yin(kin) - m * xin(kin);

            y(kout) = m * xo + b;

        elseif kin > 1

            % If sigma = 0 and xo lies beyond the final input point,
            % linearly extrapolate from the final two points.

            m = (yin(kin) - yin(kin - 1)) / ...
                (xin(kin) - xin(kin - 1));

            b = yin(kin) - m * xin(kin);

            y(kout) = m * xo + b;
        end
    end
end


function out = convolveGaussianPoint(xin, yin, k, n, xo, limit, sigma)

    SQRT2 = 1.41421356237309504880;
    SQRT2PI = 2.50662827463100050241;

    % Precompute quantities that are constant for this output point.
    %
    % These would otherwise be recalculated for every input interval.
    invSqrt2Sigma = 1.0 / (SQRT2 * sigma);
    sigmaOverSqrt2Pi = sigma / SQRT2PI;

    % 2*sigma^2, used in the Gaussian exponent.
    twoSigmaSq = 2.0 * sigma * sigma;


    % ---------------------------------------------------------------
    % Initialise at the first input point.
    % ---------------------------------------------------------------

    % Distance from the input point to the output position.
    z = xo - xin(k);

    % Unnormalised Gaussian at this point:
    %
    %   G = exp(-(xo-x)^2/(2*sigma^2))
    %
    Glo = exp(-z * z / twoSigmaSq);

    % The integral of the Gaussian is expressed in terms of erf:
    %
    %   erf(-(xo-x)/(sqrt(2)*sigma))
    %
    erfmin = erf(-z * invSqrt2Sigma);
    erflo = erfmin;

    % Accumulate the convolution integral here.
    y = 0;


    % ---------------------------------------------------------------
    % Integrate over the piecewise-linear input function.
    % ---------------------------------------------------------------

    while k < n

        k = k + 1;

        % Ignore duplicate input points.
        if xin(k) ~= xin(k - 1)

            % Distance from the new input point to xo.
            zhi = xo - xin(k);

            % Dimensionless distance in units of sigma:
            %
            %   u = -(xo-x)/(sqrt(2)*sigma)
            %
            u = -zhi * invSqrt2Sigma;

            % Unnormalised Gaussian at the new endpoint.
            Ghi = exp(-u * u);

            % Error-function value at the new endpoint.
            erfhi = erf(u);


            % -------------------------------------------------------
            % Linear interpolation between the two input points.
            %
            %   yin(x) = m*x + b
            % -------------------------------------------------------

            m = (yin(k) - yin(k - 1)) / ...
                (xin(k) - xin(k - 1));

            b = yin(k) - m * xin(k);


            % -------------------------------------------------------
            % Analytically integrate
            %
            %   (m*x + b) * Gaussian(x)
            %
            % over this input interval.
            %
            % The first term comes from the integral of the
            % Gaussian and therefore contains erf().
            %
            % The second term comes from the integral of
            % (x-xo)*Gaussian and therefore contains Ghi-Glo.
            % -------------------------------------------------------

            y = y ...
                + 0.5 * (m * xo + b) * (erfhi - erflo) ...
                - sigmaOverSqrt2Pi * m * (Ghi - Glo);


            % Current endpoint becomes the lower endpoint for
            % the next interval.
            Glo = Ghi;
            erflo = erfhi;


            % Stop once we have reached the upper Gaussian limit.
            %
            % limit ~= 3.717*sigma, corresponding to a Gaussian
            % amplitude of approximately 0.1% of its peak.
            if xin(k) >= xo + limit
                break
            end
        end
    end


    % ---------------------------------------------------------------
    % Normalisation
    % ---------------------------------------------------------------
    %
    % The Gaussian has been truncated to approximately +/-3.717*sigma,
    % so its integrated area is slightly less than one.
    %
    % The erf difference gives the area of this truncated Gaussian.
    % Renormalising ensures that a constant input function remains
    % constant after convolution.
    % ---------------------------------------------------------------

    out = 2 * y / (erflo - erfmin);

end
