drop table edges;
drop table vertex_annotations;

create table edges (prefix bigint unsigned not null,suffix bigint unsigned not null);
create index i1 on edges(prefix);
create index i2 on edges(suffix);
create table vertex_annotations(vertex bigint unsigned not null,readNumber int unsigned not null,readStrand char not null,
	readPosition smallint unsigned not null);

create index i3 on vertex_annotations(vertex);

