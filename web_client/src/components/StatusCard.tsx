import { Paper, Text, Group, Button, Loader, Stack, Tooltip } from '@mantine/core';
import { IconCheck, IconAlertTriangle, IconX } from '@tabler/icons-react';
import { JobStage, JobStats, JobStatusResponse } from '../api';

const STAGE_LABELS: Record<JobStage, string> = {
    loading: 'Loading building blocks',
    fragmenting: 'Fragmenting molecule',
    enumerating: 'Enumerating',
    profiling: 'Profiling properties',
};

interface StatusCardProps {
    status?: JobStatusResponse['status'];
    stage?: JobStage;
    stats?: JobStats;
    error?: string;
    canCancel: boolean;
    isCancelling: boolean;
    onCancel: () => void;
}

function runningLabel(status?: JobStatusResponse['status'], stage?: JobStage) {
    if (status === 'PENDING') return 'Queued — waiting for a worker';
    if (stage) return `${STAGE_LABELS[stage]}…`;
    return 'Running…';
}

export function StatusCard({ status, stage, stats, error, canCancel, isCancelling, onCancel }: StatusCardProps) {
    if (!status) return null;

    const isRunning = status !== 'SUCCESS' && status !== 'FAILURE' && status !== 'CANCELLED';

    return (
        <Paper p="md" withBorder shadow="sm">
            <Stack gap="xs">
                {isRunning && (
                    <Group justify="space-between" wrap="nowrap">
                        <Group gap="xs" wrap="nowrap">
                            <Loader size="xs" />
                            <Text size="sm">{runningLabel(status, stage)}</Text>
                        </Group>
                        <Tooltip label={canCancel ? 'Cancel job' : 'Cancel not available in local mode'} withArrow>
                            <Button
                                size="xs"
                                color="red"
                                variant="light"
                                leftSection={<IconX size={14} />}
                                disabled={!canCancel}
                                loading={isCancelling}
                                onClick={onCancel}
                            >
                                Cancel
                            </Button>
                        </Tooltip>
                    </Group>
                )}

                {status === 'SUCCESS' && (
                    <Group gap="xs" wrap="nowrap">
                        <IconCheck size={16} color="var(--mantine-color-green-6)" />
                        <Text size="sm">
                            {stats && stats.n_molecules > 0 ? (
                                <>
                                    <strong>{stats.n_molecules.toLocaleString()}</strong> molecules
                                    {stats.seconds !== null && stats.seconds !== undefined && ` in ${stats.seconds}s`}
                                </>
                            ) : (
                                'No new molecules — only the input was returned. Try adjusting the parameters.'
                            )}
                        </Text>
                    </Group>
                )}

                {status === 'FAILURE' && (
                    <Group gap="xs" wrap="nowrap" align="flex-start">
                        <IconAlertTriangle size={16} color="var(--mantine-color-red-6)" style={{ flexShrink: 0, marginTop: 2 }} />
                        <Text size="sm" style={{ wordBreak: 'break-word' }}>
                            Enumeration failed: {error || 'unknown error'}
                        </Text>
                    </Group>
                )}

                {status === 'CANCELLED' && (
                    <Group gap="xs" wrap="nowrap">
                        <IconX size={16} color="var(--mantine-color-orange-6)" />
                        <Text size="sm">Job cancelled.</Text>
                    </Group>
                )}
            </Stack>
        </Paper>
    );
}
