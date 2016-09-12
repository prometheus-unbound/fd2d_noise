
%==========================================================================
% plot forward wavefield
%==========================================================================

subplot(2, 2, 1)
max_M_tn = max(max_M_tn, max(max(abs(real(M_tn + eps)))));

if (t(n) >= 0)
    pcolor(ax1, X / 1000, Z / 1000, real(M_tn)' / max_M_tn);
else
    pcolor(ax1, X / 1000, Z / 1000, 0 * u' + eps);
end

caxis(ax1, [-0.15, 0.15])
shading interp


%- plot array ---------------------------------------------
plot(ax1, ref_station(:, 1) / 1000, ref_station(:, 2) / 1000, ...
    'kx', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)
plot(ax1, rec(:, 1) / 1000, rec(:, 2) / 1000, ...
    'kd', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)


%- plot absorbing boundaries ------------------------------
if(absorb_bottom); plot(ax1, [absorb_left*width, Lx - absorb_right*width] / 1000, [width, width] / 1000, 'k--'); end
if(absorb_top); plot(ax1, [absorb_left*width, Lx - absorb_right*width] / 1000, [Lz - width, Lz - width] / 1000, 'k--'); end
if(absorb_left); plot(ax1, [width, width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end
if(absorb_right); plot(ax1, [Lx - width, Lx - width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end


%- set labels ---------------------------------------------
xlabels = get(ax1, 'XTick');
ylabels = get(ax1, 'YTick');
set(ax1, 'XTick', [xlabels(1), xlabels(ceil(length(xlabels) / 2)), xlabels(end)])
set(ax1, 'YTick', [ylabels(1), ylabels(ceil(length(ylabels) / 2)), ylabels(end)])
xlim(ax1, [0, Lx / 1000])
ylim(ax1, [0, Lz / 1000])



%==========================================================================
% plot adjoint wavefield
%==========================================================================

subplot(2, 2, 2)
max_u = max(max_u, max(max(abs(u + eps))));

if (t(n) >= 0)
    pcolor(ax2, X / 1000, Z / 1000, u' / max_u);
elseif any(max(adjstf(:, 1:n), [], 2) > 0.07 * max(adjstf(:, 1:n_zero), [], 2))
    pcolor(ax2, X / 1000, Z / 1000, u' / max_u);
else
    pcolor(ax2, X / 1000, Z / 1000, 0 * u' + eps);
end

caxis(ax2, [-0.25, 0.25])
shading interp


%- plot array ---------------------------------------------
plot(ax2, ref_station(:, 1) / 1000, ref_station(:, 2) / 1000, ...
    'kx', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)
plot(ax2, rec(:, 1) / 1000, rec(:, 2) / 1000, ...
    'kd', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)


%- plot absorbing boundaries ------------------------------
if(absorb_bottom); plot(ax2, [absorb_left*width, Lx - absorb_right*width] / 1000, [width, width] / 1000, 'k--'); end
if(absorb_top); plot(ax2, [absorb_left*width, Lx - absorb_right*width] / 1000, [Lz - width, Lz - width] / 1000, 'k--'); end
if(absorb_left); plot(ax2, [width, width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end
if(absorb_right); plot(ax2, [Lx - width, Lx - width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end


%- set labels ---------------------------------------------
xlabels = get(ax2, 'XTick');
set(ax2, 'XTick', [xlabels(1), xlabels(ceil(length(xlabels) / 2)), xlabels(end)])
xlim(ax2, [0, Lx / 1000])
ylim(ax2, [0, Lz / 1000])



%==========================================================================
% plot kernel
%==========================================================================

subplot(2, 2, 3:4)
if (t(n) >= 0)
    pcolor(ax3, X / 1000, Z / 1000, K_source'/max(max(abs(K_source))))
    caxis(ax3, [-1, 1]);
else
    pcolor(ax3, X / 1000, Z / 1000, 0 * K_source' + eps);
end

shading interp


%- plot array ---------------------------------------------
plot(ax3, ref_station(:, 1) / 1000, ref_station(:, 2) / 1000, ...
    'kx', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)
plot(ax3, rec(:, 1) / 1000, rec(:, 2) / 1000, ...
    'kd', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)


%- plot absorbing boundaries ------------------------------
if(absorb_bottom); plot(ax3, [absorb_left*width, Lx - absorb_right*width] / 1000, [width, width] / 1000, 'k--'); end
if(absorb_top); plot(ax3, [absorb_left*width, Lx - absorb_right*width] / 1000, [Lz - width, Lz - width] / 1000, 'k--'); end
if(absorb_left); plot(ax3, [width, width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end
if(absorb_right); plot(ax3, [Lx - width, Lx - width] / 1000, [absorb_bottom*width, Lz - absorb_top*width] / 1000, 'k--'); end


%- set labels ---------------------------------------------
xlabels = get(ax3, 'XTick');
ylabels = get(ax3, 'YTick');
set(ax3, 'XTick', [xlabels(1), xlabels(ceil(length(xlabels) / 2)), xlabels(end)])
set(ax3, 'YTick', [ylabels(1), ylabels(ceil(length(ylabels) / 2)), ylabels(end)])
xlim(ax3, [0, Lx / 1000])
ylim(ax3, [0, Lz / 1000])

