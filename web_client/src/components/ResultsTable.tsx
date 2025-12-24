import { Table, Badge, ScrollArea, Text, Group, Image, Loader, Stack, Alert, ColorSwatch, Anchor } from '@mantine/core';
import { IconInfoCircle } from '@tabler/icons-react';
import { useState, useEffect, useMemo } from 'react';
import { EnumerationResult, getMolResult } from '../api';

interface ResultsTableProps {
    results: EnumerationResult[];
    jobId?: string | null;
}

const STOPLIGHT_COLORS: Record<string, string> = {
    green: 'rgba(72, 173, 85, 0.06)',
    yellow: 'rgba(235, 235, 52, 0.06)',
    red: 'rgba(235, 64, 52, 0.06)',
};

function MoleculeImage({ smiles, bbs, alpha, bgColor }: { smiles: string, bbs: string[], alpha: number, bgColor: string }) {
    const [svg, setSvg] = useState<string | null>(null);
    const [loading, setLoading] = useState(true);

    useEffect(() => {
        let mounted = true;
        setLoading(true);
        getMolResult(smiles, bbs, alpha, bgColor).then(data => {
            if (mounted) {
                setSvg(data);
                setLoading(false);
            }
        }).catch(() => {
            if (mounted) setLoading(false);
        });
        return () => { mounted = false; };
    }, [smiles, JSON.stringify(bbs), alpha, bgColor]);

    if (loading) return <Loader size="xs" />;
    if (!svg) return <Text size="xs" c="dimmed">Error</Text>;

    return <Image src={svg} w={200} />;
}

function ResultRow({ row, bbKeys }: { row: EnumerationResult, bbKeys: string[] }) {
    const [showImage, setShowImage] = useState(false);
    
    // Collect BB smiles for highlighting
    const bbs = useMemo(() => bbKeys.map(k => row[k]).filter(Boolean), [row, bbKeys]);
    
    const colorKey = row.stoplight_color ? row.stoplight_color.toLowerCase() : '';
    const bgColor = STOPLIGHT_COLORS[colorKey] || 'rgba(255, 255, 255, 1.0)';

    return (
        <Table.Tr bg={bgColor}>
            <Table.Td 
                onClick={() => setShowImage(!showImage)} 
                style={{ cursor: 'pointer' }}
            >
                {showImage ? (
                    <MoleculeImage smiles={row.Product} bbs={bbs} alpha={0.8} bgColor={bgColor} />
                ) : (
                    <>
                        <Text size="xs" c="dimmed" fs="italic">Click to view structure</Text>
                        <Text size="xs" mt={4} style={{ fontFamily: 'monospace' }}>{row.Product}</Text>
                    </>
                )}
            </Table.Td>
            <Table.Td style={{ minWidth: 80 }} ta="center">{row.Similarity_to_query?.toFixed(2)}</Table.Td>
            <Table.Td style={{ minWidth: 80 }} ta="center">{row.QED?.toFixed(2)}</Table.Td>
            <Table.Td style={{ whiteSpace: 'normal', minWidth: '180px' }}>
                {row.Reaction_name && (
                    <Badge 
                        variant="light" 
                        color="blue" 
                        style={{ 
                            textTransform: 'none', 
                            height: 'auto', 
                            padding: '4px 8px', 
                            whiteSpace: 'normal',
                            display: 'block',
                            lineHeight: '1.4',
                            textAlign: 'left'
                        }}
                    >
                        {row.Reaction_name}
                    </Badge>
                )}
            </Table.Td>
            {bbKeys.map(key => {
                const urlKey = key.replace('BB', 'URL');
                const url = row[urlKey] || (key === 'BB' ? row['URL'] : undefined);
                return (
                    <Table.Td key={key} style={{ fontFamily: 'monospace', fontSize: '12px', maxWidth: '150px', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }} title={row[key]}>
                        {url ? (
                            <Anchor href={url} target="_blank" size="xs" inherit style={{ textDecoration: 'underline' }}>
                                {row[key]}
                            </Anchor>
                        ) : (
                            row[key]
                        )}
                    </Table.Td>
                );
            })}
        </Table.Tr>
    );
}

export function ResultsTable({ results }: ResultsTableProps) {
    if (!results || results.length === 0) {
        return <Text ta="center" c="dimmed" my="xl">No results to display</Text>;
    }

    const maxRows = parseInt(import.meta.env.VITE_MAX_DISPLAY_ROWS || '250');
    const displayResults = results.slice(0, maxRows);

    const fixedKeys = ['Product', 'Similarity_to_query', 'Reaction_name', 'stoplight_color', 'QED'];
    
    const bbKeys = useMemo(() => {
        // Scan all results to find all unique BB keys
        const allKeys = new Set<string>();
        displayResults.forEach(row => {
            Object.keys(row).forEach(k => {
                if (k.startsWith('BB') && !fixedKeys.includes(k)) {
                    allKeys.add(k);
                }
            });
        });

        return Array.from(allKeys).sort((a, b) => {
            const numA = parseInt(a.replace('BB', '')) || 0;
            const numB = parseInt(b.replace('BB', '')) || 0;
            return numA - numB;
        });
    }, [displayResults]);

    return (
        <Stack gap="sm">
            <Alert icon={<IconInfoCircle size={16} />} title="Results Legend" color="gray" variant="light">
                 <Group gap="md">
                    <Group gap="xs"><ColorSwatch color="rgba(72, 173, 85, 0.6)" size={14} /> <Text size="sm">High Drug-likeness (Green)</Text></Group>
                    <Group gap="xs"><ColorSwatch color="rgba(235, 235, 52, 0.6)" size={14} /> <Text size="sm">Medium Drug-likeness (Yellow)</Text></Group>
                    <Group gap="xs"><ColorSwatch color="rgba(235, 64, 52, 0.6)" size={14} /> <Text size="sm">Low Drug-likeness (Red)</Text></Group>
                 </Group>
                 <Text size="xs" mt="xs">
                    The table below displays the first {maxRows} enumeration results at max. Columns are also truncated for simplicity.
                    To view the complete data, download the full results. Rows are colored based on{' '}
                    <Anchor href="https://doi.org/10.1021/acs.jcim.4c00412" target="_blank" inherit>
                        Stoplight Scores
                    </Anchor>
                    , which is a drug-likeness indicator. Click on a row to visualize the molecule. The building blocks are highlighted 
                    in the structure but may not be accurate always. You can also click on building block IDs to visit their source page 
                    (if available).
                 </Text>
            </Alert>
            
            <ScrollArea h={600}>
                 <Table striped highlightOnHover verticalSpacing="sm">
                    <Table.Thead>
                        <Table.Tr>
                            <Table.Th w={220}>Product</Table.Th>
                            <Table.Th w={80} ta="center">Sim.</Table.Th>
                            <Table.Th w={80} ta="center">QED</Table.Th>
                            <Table.Th w={200}>Reaction Path</Table.Th>
                            {bbKeys.map(key => <Table.Th key={key}>{key}</Table.Th>)}
                        </Table.Tr>
                    </Table.Thead>
                    <Table.Tbody>
                        {displayResults.map((row, index) => (
                            <ResultRow key={index} row={row} bbKeys={bbKeys} />
                        ))}
                    </Table.Tbody>
                </Table>
            </ScrollArea>
        </Stack>
    );
}
