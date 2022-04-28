library(ggplot2)

## Fingerprint vector comparison

count = read.csv('morgan_count.csv', header = FALSE)
sim = read.csv('morgan_sim.csv', header = FALSE)
bit = read.csv('morgan_bit.csv', header = FALSE)

count =cbind('Count', count)
colnames(count) = c('Vector','Reaction_pair','Score')

sim = cbind('Count simulation', sim)
colnames(sim) = c('Vector','Reaction_pair','Score')

bit = cbind('No count simulation', bit)
colnames(bit) = c('Vector','Reaction_pair','Score')

bit_vector_types = rbind(bit, sim)
scatter = merge(bit_vector_types, count, by = 'Reaction_pair')
colnames(scatter) = c('Reaction_pair', 'Bit_vect_type', 'Bit_score', 'Count_vect', 'Count_score')

p =ggplot(scatter, aes(x=Count_score, y=Bit_score, color=Bit_vect_type))
p+geom_point(size = 0.6)+  theme(legend.key.size = unit(2, 'cm'),legend.text = element_text(size=12), legend.title = element_text(size=0)) + 
  xlab('Count vector similarity score') + ylab('Bit vector similarity score')+ guides(colour = guide_legend(override.aes = list(size=1.5)))+ 
  theme(legend.position = 'top')
cor(scatter$Count_score[scatter$Bit_vect_type=='Count simulation'], scatter$Bit_score[scatter$Bit_vect_type=='Count simulation'])
cor(scatter$Count_score[scatter$Bit_vect_type=='No count simulation'], scatter$Bit_score[scatter$Bit_vect_type=='No count simulation'])

## Dumbell plot reaction centre methods

mapping_sims = read.csv('similarity_ap_mapping.csv', header = FALSE) # or m_mapping/r_mapping
mapping_sims = cbind(mapping_sims, 'Atom_mapping')
colnames(mapping_sims) = c('Pair', 'Score', 'Method')

subtraction_sims = read.csv('similarity_ap_subtraction.csv', header = FALSE) #or m_subtraction/r_subtraction
subtraction_sims = cbind(subtraction_sims, 'Subtraction')
colnames(subtraction_sims) = c('Pair', "Score", 'Method')

sims = rbind(mapping_sims, subtraction_sims)
ggplot(data = sims,aes(x= Score, y= reorder(Pair, Score))) +
geom_line(aes(group = Pair))+
geom_point(aes(color=Method), size=4) +labs(x = 'Similarity score', y = 'Reaction pair') +
scale_color_manual(labels = c("Atom-atom mapping", "Fingerprint subtraction"), values = c("#F8766D", "#00BFC4")) +
xlim(0.16, 1.0)

## Accuracy at different transformation-level similarity weightings

range = 0:100
reaction_numbers = 0:19
same_enz = read.csv('accuracy_subtraction_m.csv')  # or subtraction_r/subtraction_ap/mapping_m/mapping_r etc.
percent_correct = data.frame(weighting = numeric(),
                             accuracy = numeric())

for (number in range)
{
  weighting = number/100
  same_enz$overall_score = ((same_enz$transformation_score*weighting) + (same_enz$whole_mol_score*(1-weighting)))
  correctly_identified = 0
  for (reaction in reaction_numbers)
  {
    if (same_enz$overall_score[same_enz$test_reaction==same_enz$correct_reaction & same_enz$reaction==reaction]==max(same_enz$overall_score[same_enz$reaction==reaction]))
    {
      correctly_identified = correctly_identified +1
    }
  }
  percent_correct = rbind(percent_correct,c(weighting, correctly_identified/20))
}

colnames(percent_correct) = c('Weighting','Percent_correct')
q = ggplot(percent_correct, aes(x=Weighting*100, y=Percent_correct*100)) + geom_line() + labs(y = '% of reactions successfully identified', x = 'Transformation-level similarity weighting (%)')
q

#Fragment size plot

subtraction_frags = read.csv('subtraction_fragment_size.csv')

ggplot(data = subtraction_frags, aes(x=Size, y=Rank, group = Fingerprint)) + geom_line(aes(color = Fingerprint, linetype = Fingerprint)) + ylim(5,1)+
  geom_point(aes(color = Fingerprint)) +labs(x = 'Fragment Size', y = 'Rank')+
  scale_color_manual(labels = c("Atom Pair", "Morgan", 'Topological'), values = c("#F8766D", 'chartreuse3',"#00BFC4")) 

mapping_frags = read.csv('Mapping_fragment_size.csv')
ggplot(data = mapping_frags, aes(x=Size, y=Rank, group = Fingerprint)) + geom_line(aes(color = Fingerprint, linetype = Fingerprint)) + ylim(30,1)+
  geom_point(aes(color = Fingerprint)) +labs(x = 'Distance from changing atoms (number of bonds)', y = 'Rank')+
  scale_color_manual(labels = c("Atom Pair", "Morgan", 'Topological'), values = c("#F8766D", 'chartreuse3',"#00BFC4")) 
    