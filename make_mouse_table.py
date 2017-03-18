from basic import *
import html_colors
import util

with Parser(locals()) as p:
#    p.str('args').unspecified_default().multiple().required()
    p.str('clones_file').required()
    p.str('outfile_prefix')
    p.flag('horizontal_lines')
    p.flag('show')
    p.flag('include_counts_in_mouse_labels')

if not outfile_prefix:
    outfile_prefix = clones_file[:-4]

import matplotlib
if not show: matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import numpy as np




# all_tcrs = {}
# all_epitopes = []

# infields = []
# for line in open( clones_file,'r'):
#     if not infields:
#         if line[0] == '#':
#             infields = line[1:-1].split('\t')
#         else:
#             infields = line[:-1].split('\t')
#         continue
#     assert infields

#     l = parse_tsv_line( line[:-1], infields )

#     mouse = l['mouse']
#     epitope = l['epitope']
#     clone_size = int(l['clone_size'])

#     if mouse not in all_tcrs:
#         all_tcrs[mouse] = {}
#     if epitope not in all_tcrs[mouse]:
#         all_tcrs[mouse][epitope] = []
#         if epitope not in all_epitopes:
#             all_epitopes.append( epitope )
#     all_tcrs[mouse][epitope].append( clone_size ) ## just store the clone sizes


all_tcrs = parse_tsv_file( clones_file, ['subject','epitope'], ['clone_size'], False )
all_epitopes = list( reduce( set.union, ( set( x.keys() ) for x in all_tcrs.values() ) ) )
all_epitopes.sort()

all_mice= all_tcrs.keys()[:]
all_mice.sort()

counts = {}
for e in all_epitopes: counts[e] = [0,0]
for m in all_mice: counts[m] = [0,0]

for mouse in all_tcrs:
    for epitope in all_tcrs[mouse]:
        clone_sizes = [int(x[0]) for x in all_tcrs[mouse][epitope]]
        total_reads = sum(clone_sizes)
        for k in [mouse,epitope]:
            counts[k][0] += len(clone_sizes)
            counts[k][1] += total_reads

mouse_labels = {}
for mouse in all_mice:
    if include_counts_in_mouse_labels:
        mouse_labels[mouse] = '{} ({};{})'.format( mouse, counts[mouse][0], counts[mouse][1] )
    else:
        mouse_labels[mouse] = mouse

epitope_labels = {}
for epitope in all_epitopes:
    epitope_labels[epitope] = '{}  ({};{})'.format( epitope, counts[epitope][0], counts[epitope][1] )

nrows = len( all_mice )
ncols = len( all_epitopes )

preferred_plot_width = 12.0
preferred_plot_height = 12.0

preferred_cell_size = max( 0.5, min( preferred_plot_height/nrows, preferred_plot_width/ncols ) )

plot_width = ncols * preferred_cell_size
plot_height = nrows * preferred_cell_size

fontsize_small = 8.
fontsize_medium = 10.
fontsize_names = 12.

for repeat in range(3):
    if plot_width <= 1.2 * preferred_plot_width and plot_height <= 1.2 * preferred_plot_height: break

    if plot_width / preferred_plot_width > plot_height / preferred_plot_height: ## too wide
        plot_width *= 0.75
        plot_height *= 0.9
        fontsize_small *= 0.9
        fontsize_medium *= 0.9

    else: ## too tall
        plot_height *= 0.75
        plot_width *= 0.9
        fontsize_small *= 0.9
        fontsize_medium *= 0.9

fontsize_small = max(5,int(floor(0.5+fontsize_small)))
fontsize_medium = max(6,int(floor(0.5+fontsize_medium)))


fudge = 1.2
bottom_spacer = 0.3 # inches
left_margin_inches   = fudge * max( ( len(mouse_labels[x]) for x in all_mice ) ) * 0.6 * fontsize_names / 72.0
bottom_margin_inches = fudge * max( ( len(epitope_labels[x]) for x in all_epitopes ) ) * 0.75 * fontsize_names / 72.0 + bottom_spacer

top_margin_inches = 0.25
right_margin_inches = 0.25

fig_width  =   left_margin_inches + plot_width  + right_margin_inches
fig_height = bottom_margin_inches + plot_height +   top_margin_inches

top_margin    = float( bottom_margin_inches + plot_height ) / fig_height
bottom_margin = float( bottom_margin_inches ) / fig_height
left_margin   = float( left_margin_inches ) / fig_width
right_margin  = float( left_margin_inches + plot_width ) / fig_width


print 'fig_width: {:.1f} fig_height: {:.1f}'.format(fig_width,fig_height)

fig = plt.figure(1,figsize=(fig_width,fig_height))

#fig = plt.figure(1,figsize=(23,8))

#fig1.add_line(Line2D([0.5,0.5], [0,1], linewidth=2, color='blue'))
#ax = fig.add_axes( [ left_margin, bottom_margin, right_margin,top_margin ] )
#ax.grid(True)




plotno=0
for mouse in all_mice:
    for epitope in all_epitopes:
        plotno += 1
        if epitope not in all_tcrs[mouse]:
            continue

        plt.subplot( nrows, ncols, plotno )

        clone_sizes = [int(x[0]) for x in all_tcrs[mouse][epitope]]
        clone_sizes.sort()
        clone_sizes.reverse()

        colors = html_colors.get_rank_colors_no_lights(len(clone_sizes))


        wedges, texts = plt.pie( clone_sizes )
        for ii,w in enumerate(wedges):
            w.set_edgecolor('none')
            w.set_facecolor(colors[ii])

        topsize = clone_sizes[0]
        total_size = sum(clone_sizes)

        ## show the size of the largest wedge?
        if len(wedges)>1:
            w = wedges[0]
            #print w.center, w.r, w.theta1, w.theta2
            ## show the size at radius distance in middle of edge
            angle_degrees = w.theta2*0.5
            if 65<=angle_degrees<=115: angle_degrees = 65. if angle_degrees < 90. else 115.
            x=1.1*w.r*math.cos( math.pi * angle_degrees / 180.0 )
            y=1.1*w.r*math.sin( math.pi * angle_degrees / 180.0 )
            thresh = 0.3*w.r
            ha = 'left'   if x>thresh else ( 'center' if x>-thresh else 'right' )
            va = 'bottom' if y>thresh else ( 'center' if y>-thresh else 'top' )
            plt.text(x,y,`topsize`,fontdict={'fontsize':fontsize_small},color='r',
                     horizontalalignment=ha,verticalalignment=va)


        ## show the total number of reads
        radius = wedges[0].r
        plt.text(0,-1.1*radius,`total_size`,fontdict={'fontsize':fontsize_medium},
                 horizontalalignment='center',verticalalignment='top' )

        #t = plt.title(`sum(clone_sizes)`,fontdict={'fontsize':8})
        if False:
            if epitope==all_epitopes[0]:
                plt.title(mouse)
            elif mouse==all_mice[0]:
                plt.title(epitope)
        #break
    #break


#plt.hlines(0.5,0.0,1.0)
#plt.vlines(0.5,0.0,1.0)


epsilon = 0.0
plt.subplots_adjust(
    left=left_margin+epsilon,
    right=right_margin-epsilon,
    bottom=bottom_margin+epsilon,
    top=top_margin-epsilon
)

ywidth = (top_margin-bottom_margin) / ( len(all_mice) )
xwidth = (right_margin-left_margin) / ( len(all_epitopes) )

#ystep = (top_margin-bottom_margin) / ( len(all_epitopes)-1 )
lines = []


# if horizontal_lines:
#     for ii in range(len(all_epitopes)):
#     #for ii in range(len(all_epitopes)+1):
#         y = bottom_margin + 1.02 * ii * ywidth
#         lines.append( matplotlib.lines.Line2D( [0,1], [y,y],
#                                                transform=fig.transFigure, figure=fig, c='k' ) )

if False:
    for ii in range(len(all_mice)+1):
        x = left_margin + ii*xwidth
        lines.append( matplotlib.lines.Line2D( [x,x], [0,1],
                                               transform=fig.transFigure, figure=fig, c='k' ) )

    fig.lines.extend(lines)

for ii,mouse in enumerate( all_mice ):
    plt.figtext( left_margin-0.005, top_margin - 3*ywidth/5 - ii * ywidth, mouse_labels[mouse], ha='right', va='center',
                 fontdict={'fontsize':fontsize_names})
    #plt.figtext( right_margin+0.005, top_margin - 3*ywidth/5 - ii * ywidth, epitope,ha='left')

#xstep = (right_margin-left_margin) / ( len(all_mice)-1 )

for ii,epitope in enumerate( all_epitopes ):
    #name = mouse[:]
    # if name[0] == 'd' and 'Mouse' in name:
    #     name = name.replace('Mouse','_')
    plt.figtext(left_margin + xwidth/2 + ii * xwidth, bottom_margin - (bottom_spacer)/fig_height,
                epitope_labels[epitope],
                rotation='vertical', ha='center', va='top',
                fontdict={'fontsize':fontsize_names})

    #plt.figtext(left_margin + xwidth/2 + ii * xwidth, 0.98, epitope, ha='center', va='top' )



pngfile = outfile_prefix+'_subject_table.png'
print 'making:',pngfile
plt.savefig(pngfile)

util.readme(pngfile,"""This subject-table plot shows all the successfully parsed, paired reads, split by mouse/subject (the rows)
and epitope (the columns, labeled at the bottom). The epitope column labels include in parentheses the number of clones followed by
the total number of TCRs. Each pie shows the paired reads for a single mouse/epitope combination, with each wedge corresponding to
a clone. The size of the top clone is shown in red near the red wedge, and the total number of reads is shown below the pie in black.
""")


if show:
    plt.show()


